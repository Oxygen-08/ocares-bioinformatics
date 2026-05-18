#!/usr/bin/env python3
"""
Run NUCmer on 5 strains not yet in tiered_blocks.tsv:
  - RS218 and APEC_O1 (newly downloaded)
  - HUSEC2011, NRG857C, TW10598 (already in manifest but never aligned)
Then APPEND their alignment blocks to the existing tiered_blocks.tsv.

Run from repo root:
  /opt/anaconda3/envs/fp_pipeline/bin/python scripts/analysis/run_new_pathogens_nucmer.py
"""

import csv
import logging
import subprocess
from dataclasses import asdict, dataclass
from pathlib import Path

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s  %(levelname)s  %(message)s",
    datefmt="%H:%M:%S",
)
log = logging.getLogger(__name__)

NUCMER       = Path("/opt/anaconda3/envs/fp_pipeline/bin/nucmer")
DELTA_FILTER = Path("/opt/anaconda3/envs/fp_pipeline/bin/delta-filter")
SHOW_COORDS  = Path("/opt/anaconda3/envs/fp_pipeline/bin/show-coords")

REPO_ROOT   = Path(__file__).parents[2]
GENOME_DIR  = REPO_ROOT / "data" / "genomes"
RESULTS_DIR = REPO_ROOT / "data" / "results" / "nucmer"
BLOCKS_TSV  = REPO_ROOT / "data" / "results" / "tiered_blocks.tsv"
REFERENCE   = GENOME_DIR / "O157H7_Sakai" / "O157H7_Sakai.fna"

MANIFEST_PREFIX = "/Users/user/bioinformatics"
ACTUAL_PREFIX   = str(REPO_ROOT)

TARGET_LABELS = [
    "HUSEC2011", "NRG857C", "TW10598",
    "RS218", "APEC_O1",
]


@dataclass
class AlignBlock:
    label:      str
    pathotype:  str
    ref_start:  int
    ref_end:    int
    qry_start:  int
    qry_end:    int
    ref_len:    int
    qry_len:    int
    identity:   float
    ref_contig: str
    qry_contig: str
    tier:       str


def classify_tier(identity: float) -> str:
    if identity >= 95.0:
        return "CONSERVED"
    if identity >= 85.0:
        return "MODERATE"
    return "DIVERGED"


def run_nucmer(query: Path, prefix: Path) -> Path:
    delta = prefix.with_suffix(".delta")
    if delta.exists():
        log.info("  SKIP (delta exists): %s", delta.name)
        return delta
    log.info("  Running nucmer: %s", query.parent.name)
    subprocess.run(
        [str(NUCMER), "--maxgap=500", "--mincluster=65",
         "--prefix", str(prefix), str(REFERENCE), str(query)],
        check=True, capture_output=True,
    )
    return delta


def filter_delta(delta: Path) -> Path:
    filtered = delta.with_suffix(".filtered.delta")
    if filtered.exists():
        return filtered
    with open(filtered, "w") as fh:
        subprocess.run(
            [str(DELTA_FILTER), "-1", str(delta)],
            stdout=fh, check=True, capture_output=False,
        )
    return filtered


def show_coords(filtered: Path) -> Path:
    coords = filtered.with_suffix(".coords")
    if coords.exists():
        return coords
    with open(coords, "w") as fh:
        subprocess.run(
            [str(SHOW_COORDS), "-T", "-r", "-l", str(filtered)],
            stdout=fh, check=True, capture_output=False,
        )
    return coords


def parse_coords(coords_file: Path, label: str, pathotype: str) -> list[AlignBlock]:
    blocks = []
    with open(coords_file) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("[") or line.startswith("="):
                continue
            parts = line.split("\t")
            if len(parts) < 11:
                continue
            try:
                r1, r2, q1, q2, rlen, qlen, ident, ref_ctg, qry_ctg = (
                    int(parts[0]), int(parts[1]),
                    int(parts[2]), int(parts[3]),
                    int(parts[4]), int(parts[5]),
                    float(parts[6]),
                    parts[9], parts[10],
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


def load_manifest() -> dict[str, tuple[str, str]]:
    manifest = GENOME_DIR / "genome_manifest.tsv"
    result = {}
    with open(manifest) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            if row["label"] not in TARGET_LABELS:
                continue
            fasta = row["fasta_path"].replace(MANIFEST_PREFIX, ACTUAL_PREFIX)
            result[row["label"]] = (row["pathotype"], fasta)
    return result


def main() -> None:
    # Guard: refuse to append strains already in tiered_blocks.tsv
    existing_labels: set[str] = set()
    if BLOCKS_TSV.exists():
        with open(BLOCKS_TSV) as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for row in reader:
                existing_labels.add(row["label"])

    already_done = [lbl for lbl in TARGET_LABELS if lbl in existing_labels]
    to_run = [lbl for lbl in TARGET_LABELS if lbl not in existing_labels]
    if already_done:
        log.info("Already in tiered_blocks.tsv (skip): %s", already_done)
    if not to_run:
        log.info("All target strains already present — nothing to do.")
        return

    strains = load_manifest()
    missing_from_manifest = [lbl for lbl in to_run if lbl not in strains]
    if missing_from_manifest:
        log.error("Labels not found in manifest: %s", missing_from_manifest)

    new_blocks: list[AlignBlock] = []
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)

    for label in to_run:
        if label not in strains:
            continue
        pathotype, fasta_path = strains[label]
        query = Path(fasta_path)
        if not query.exists():
            log.error("FASTA not found — skipping %s: %s", label, query)
            continue

        work_dir = RESULTS_DIR / label
        work_dir.mkdir(exist_ok=True)
        prefix = work_dir / label

        try:
            delta    = run_nucmer(query, prefix)
            filtered = filter_delta(delta)
            coords   = show_coords(filtered)
            blocks   = parse_coords(coords, label, pathotype)
        except subprocess.CalledProcessError as exc:
            log.error("Pipeline failed for %s: %s", label, exc)
            continue

        n_c = sum(1 for b in blocks if b.tier == "CONSERVED")
        n_m = sum(1 for b in blocks if b.tier == "MODERATE")
        n_d = sum(1 for b in blocks if b.tier == "DIVERGED")
        log.info("  %s (%s): %d blocks — %d C / %d M / %d D",
                 label, pathotype, len(blocks), n_c, n_m, n_d)
        new_blocks.extend(blocks)

    if not new_blocks:
        log.warning("No new blocks generated — nothing to append.")
        return

    fields = list(asdict(new_blocks[0]).keys())
    write_header = not BLOCKS_TSV.exists()

    with open(BLOCKS_TSV, "a", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
        if write_header:
            writer.writeheader()
        for b in new_blocks:
            writer.writerow(asdict(b))

    log.info("Appended %d new alignment blocks to %s", len(new_blocks), BLOCKS_TSV.name)


if __name__ == "__main__":
    main()
