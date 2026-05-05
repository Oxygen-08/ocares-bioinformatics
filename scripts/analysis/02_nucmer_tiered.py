#!/usr/bin/env python3
"""
Phase 1 — Tiered comparative genomics with NUCmer (MUMmer4).

For each non-reference genome, runs:
  nucmer  → delta file (all pairwise alignments)
  delta-filter -1 → best 1-to-1 hits (removes redundant overlapping blocks)
  show-coords -T -r -l → tabular alignment coordinates with % identity

Each alignment block is then classified into one of three tiers:

  Tier I   — Conserved         : identity >= 95.0%
  Tier II  — Moderately Diverged: 85.0% <= identity < 95.0%
  Tier III — Highly Diverged   : identity < 85.0%

Tier II and III blocks in the query (O157:H7 Sakai) are candidate
discriminatory markers — genomic regions that distinguish EHEC from
non-pathogenic relatives and commensals.

Output:
  data/results/nucmer/<label>/        delta, filtered delta, coords files
  data/results/tiered_blocks.tsv      all alignment blocks with tier label
  data/results/tier_summary.tsv       per-genome tier statistics
"""

import subprocess
import sys
import logging
from pathlib import Path
from dataclasses import dataclass, asdict
import csv

BIN_DIR = Path(sys.executable).parent

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s  %(levelname)s  %(message)s",
    datefmt="%H:%M:%S",
)
log = logging.getLogger(__name__)

REPO_ROOT    = Path(__file__).parents[2]
GENOME_DIR   = REPO_ROOT / "data" / "genomes"
RESULTS_DIR  = REPO_ROOT / "data" / "results" / "nucmer"
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

REFERENCE    = GENOME_DIR / "O157H7_Sakai" / "O157H7_Sakai.fna"

TIER_CONSERVED   = (95.0, 100.0)   # >= 95%
TIER_MODERATE    = (85.0,  95.0)   # 85 – <95%
TIER_DIVERGED    = (0.0,   85.0)   # < 85%


@dataclass
class AlignBlock:
    label:        str
    pathotype:    str
    ref_start:    int
    ref_end:      int
    qry_start:    int
    qry_end:      int
    ref_len:      int
    qry_len:      int
    identity:     float
    ref_contig:   str
    qry_contig:   str
    tier:         str   # CONSERVED | MODERATE | DIVERGED


def classify_tier(identity: float) -> str:
    if identity >= TIER_CONSERVED[0]:
        return "CONSERVED"
    if identity >= TIER_MODERATE[0]:
        return "MODERATE"
    return "DIVERGED"


def run_nucmer(ref: Path, query: Path, out_prefix: Path) -> Path:
    """Run nucmer; return path to .delta file."""
    delta = out_prefix.with_suffix(".delta")
    if delta.exists():
        log.info("  SKIP nucmer (delta exists): %s", delta.name)
        return delta

    cmd = [
        str(BIN_DIR / "nucmer"),
        "--maxgap=500",
        "--mincluster=65",
        "--prefix", str(out_prefix),
        str(ref),
        str(query),
    ]
    log.info("  nucmer %s vs reference", query.parent.name)
    subprocess.run(cmd, check=True, capture_output=True)
    return delta


def filter_delta(delta: Path) -> Path:
    """Apply delta-filter -1 (best 1-to-1 global hits) to remove tandem duplicates."""
    filtered = delta.with_suffix(".filtered.delta")
    if filtered.exists():
        return filtered
    with open(filtered, "w") as fh:
        subprocess.run(
            [str(BIN_DIR / "delta-filter"), "-1", str(delta)],
            stdout=fh, check=True, capture_output=False,
        )
    return filtered


def parse_coords(coords_file: Path, label: str, pathotype: str) -> list[AlignBlock]:
    """Parse show-coords -T -r -l output into AlignBlock list."""
    blocks = []
    with open(coords_file) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("[") or line.startswith("="):
                continue
            parts = line.split("\t")
            if len(parts) < 9:
                continue
            try:
                r1, r2, q1, q2, rlen, qlen, ident, ref_ctg, qry_ctg = (
                    int(parts[0]), int(parts[1]),
                    int(parts[2]), int(parts[3]),
                    int(parts[4]), int(parts[5]),
                    float(parts[6]),
                    parts[7], parts[8],
                )
            except (ValueError, IndexError):
                continue
            blocks.append(AlignBlock(
                label=label, pathotype=pathotype,
                ref_start=r1, ref_end=r2,
                qry_start=q1, qry_end=q2,
                ref_len=rlen, qry_len=qlen,
                identity=ident,
                ref_contig=ref_ctg, qry_contig=qry_ctg,
                tier=classify_tier(ident),
            ))
    return blocks


def show_coords(filtered_delta: Path) -> Path:
    """Run show-coords to generate tabular alignment report."""
    coords = filtered_delta.with_suffix(".coords")
    if coords.exists():
        return coords
    with open(coords, "w") as fh:
        subprocess.run(
            [str(BIN_DIR / "show-coords"), "-T", "-r", "-l", str(filtered_delta)],
            stdout=fh, check=True, capture_output=False,
        )
    return coords


def load_manifest() -> list[tuple[str, str, str]]:
    """Load genome_manifest.tsv; return list of (label, pathotype, fasta_path)."""
    manifest = GENOME_DIR / "genome_manifest.tsv"
    if not manifest.exists():
        log.error("Manifest not found — run 01_download_genomes.py first")
        sys.exit(1)
    entries = []
    with open(manifest) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            if row["fasta_path"] == "FAILED":
                log.warning("SKIP %s — download failed", row["label"])
                continue
            if row["label"] == "O157H7_Sakai":
                continue  # skip self-comparison
            entries.append((row["label"], row["pathotype"], row["fasta_path"]))
    return entries


def main() -> None:
    if not REFERENCE.exists():
        log.error("Reference genome not found: %s", REFERENCE)
        sys.exit(1)

    strains = load_manifest()
    log.info("Running NUCmer on %d strains against O157:H7 Sakai", len(strains))

    all_blocks: list[AlignBlock] = []
    summary_rows: list[dict] = []

    for label, pathotype, fasta_path in strains:
        query = Path(fasta_path)
        work_dir = RESULTS_DIR / label
        work_dir.mkdir(exist_ok=True)
        prefix = work_dir / label

        try:
            delta    = run_nucmer(REFERENCE, query, prefix)
            filtered = filter_delta(delta)
            coords   = show_coords(filtered)
            blocks   = parse_coords(coords, label, pathotype)
        except subprocess.CalledProcessError as exc:
            log.error("Pipeline failed for %s: %s", label, exc)
            continue

        all_blocks.extend(blocks)

        n_conserved = sum(1 for b in blocks if b.tier == "CONSERVED")
        n_moderate  = sum(1 for b in blocks if b.tier == "MODERATE")
        n_diverged  = sum(1 for b in blocks if b.tier == "DIVERGED")
        total_ref   = sum(b.ref_len for b in blocks)

        summary_rows.append({
            "label":       label,
            "pathotype":   pathotype,
            "n_blocks":    len(blocks),
            "n_conserved": n_conserved,
            "n_moderate":  n_moderate,
            "n_diverged":  n_diverged,
            "total_aligned_bp": total_ref,
            "pct_conserved": round(100 * n_conserved / max(len(blocks), 1), 2),
        })
        log.info("  %s: %d blocks — %d conserved / %d moderate / %d diverged",
                 label, len(blocks), n_conserved, n_moderate, n_diverged)

    # Write tiered block table
    blocks_out = REPO_ROOT / "data" / "results" / "tiered_blocks.tsv"
    if all_blocks:
        fields = list(asdict(all_blocks[0]).keys())
        with open(blocks_out, "w", newline="") as fh:
            writer = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
            writer.writeheader()
            for b in all_blocks:
                writer.writerow(asdict(b))
        log.info("Tiered blocks written: %s (%d rows)", blocks_out, len(all_blocks))

    # Write summary table
    summary_out = REPO_ROOT / "data" / "results" / "tier_summary.tsv"
    if summary_rows:
        with open(summary_out, "w", newline="") as fh:
            writer = csv.DictWriter(fh, fieldnames=summary_rows[0].keys(), delimiter="\t")
            writer.writeheader()
            writer.writerows(summary_rows)
        log.info("Tier summary written: %s", summary_out)


if __name__ == "__main__":
    main()
