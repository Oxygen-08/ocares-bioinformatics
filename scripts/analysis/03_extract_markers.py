#!/usr/bin/env python3
"""
Extract Tier II + III divergent blocks from O157:H7 Sakai as candidate
pathogen-discriminatory markers, then export as multi-FASTA for BLASTn.

Logic:
  - Only blocks where tier != CONSERVED are potential discriminatory markers
  - Merge overlapping blocks within the same tier (bedtools merge equivalent,
    implemented directly to avoid subprocess overhead on small datasets)
  - Minimum marker length: 200 bp (below this, BLASTn specificity degrades)
  - Output one FASTA per tier for parallel BLASTn screening

Output:
  data/results/markers/tier_II_markers.fna   — moderate divergence candidates
  data/results/markers/tier_III_markers.fna  — high divergence candidates
  data/results/markers/all_discriminatory.fna — combined
  data/results/markers/marker_metadata.tsv
"""

import csv
import logging
from collections import defaultdict
from pathlib import Path

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s  %(levelname)s  %(message)s",
    datefmt="%H:%M:%S",
)
log = logging.getLogger(__name__)

REPO_ROOT    = Path(__file__).parents[2]
RESULTS_DIR  = REPO_ROOT / "data" / "results"
MARKERS_DIR  = RESULTS_DIR / "markers"
MARKERS_DIR.mkdir(parents=True, exist_ok=True)

REFERENCE    = REPO_ROOT / "data" / "genomes" / "O157H7_Sakai" / "O157H7_Sakai.fna"
MIN_LEN      = 200   # bp — minimum marker length for BLASTn


def load_reference() -> dict[str, str]:
    """Load reference FASTA into {contig_id: sequence} dict."""
    log.info("Loading reference genome: %s", REFERENCE.name)
    return {rec.id: str(rec.seq) for rec in SeqIO.parse(str(REFERENCE), "fasta")}


def load_tiered_blocks() -> list[dict]:
    """Load tiered_blocks.tsv produced by 02_nucmer_tiered.py."""
    blocks_file = RESULTS_DIR / "tiered_blocks.tsv"
    if not blocks_file.exists():
        raise FileNotFoundError(f"Run 02_nucmer_tiered.py first: {blocks_file}")
    rows = []
    with open(blocks_file) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            row["ref_start"] = int(row["ref_start"])
            row["ref_end"]   = int(row["ref_end"])
            row["identity"]  = float(row["identity"])
            rows.append(row)
    return rows


def merge_intervals(intervals: list[tuple[int, int]]) -> list[tuple[int, int]]:
    """Merge overlapping or adjacent intervals; return sorted merged list."""
    if not intervals:
        return []
    sorted_iv = sorted(intervals)
    merged = [sorted_iv[0]]
    for start, end in sorted_iv[1:]:
        if start <= merged[-1][1] + 1:
            merged[-1] = (merged[-1][0], max(merged[-1][1], end))
        else:
            merged.append((start, end))
    return merged


def collect_divergent_regions(blocks: list[dict]) -> dict[str, list[tuple[int, int, str]]]:
    """
    For each reference contig, collect the union of Tier II + III intervals
    across ALL strains compared.  Regions consistently diverged across multiple
    strains are the strongest candidate markers.

    Returns: {contig_id: [(start, end, tier), ...]}
    """
    # per-tier intervals per contig
    contig_intervals: dict[str, dict[str, list[tuple[int, int]]]] = defaultdict(
        lambda: {"MODERATE": [], "DIVERGED": []}
    )

    for block in blocks:
        if block["tier"] == "CONSERVED":
            continue
        contig = block["ref_contig"]
        tier   = block["tier"]
        contig_intervals[contig][tier].append((block["ref_start"], block["ref_end"]))

    result: dict[str, list[tuple[int, int, str]]] = defaultdict(list)
    for contig, tiers in contig_intervals.items():
        for tier_label, ivs in tiers.items():
            for start, end in merge_intervals(ivs):
                if end - start + 1 >= MIN_LEN:
                    result[contig].append((start, end, tier_label))

    return result


def write_marker_fasta(
    regions: dict[str, list[tuple[int, int, str]]],
    ref_seqs: dict[str, str],
    tier_filter: str | None = None,
    out_path: Path | None = None,
) -> list[dict]:
    """Extract sequences and write FASTA; return metadata rows."""
    records = []
    metadata = []
    marker_id = 0

    for contig, intervals in regions.items():
        if contig not in ref_seqs:
            log.warning("Contig %s not found in reference — skipping", contig)
            continue
        seq_str = ref_seqs[contig]
        for start, end, tier in sorted(intervals):
            if tier_filter and tier != tier_filter:
                continue
            length = end - start + 1
            marker_id += 1
            mid   = f"MARKER_{marker_id:04d}"
            seq   = seq_str[start - 1:end]   # coords are 1-based inclusive
            header = f"{mid}|{contig}|{start}-{end}|{tier}|len={length}"
            records.append(SeqRecord(Seq(seq), id=header, description=""))
            metadata.append({
                "marker_id": mid,
                "contig":    contig,
                "start":     start,
                "end":       end,
                "length":    length,
                "tier":      tier,
            })

    if out_path and records:
        with open(out_path, "w") as fh:
            SeqIO.write(records, fh, "fasta")
        log.info("Written %d markers → %s", len(records), out_path.name)

    return metadata


def main() -> None:
    ref_seqs = load_reference()
    blocks   = load_tiered_blocks()

    log.info("Loaded %d alignment blocks", len(blocks))
    regions = collect_divergent_regions(blocks)
    total_regions = sum(len(v) for v in regions.values())
    log.info("Merged divergent regions: %d across %d contigs", total_regions, len(regions))

    # Tier II markers
    meta_ii  = write_marker_fasta(regions, ref_seqs,
                                  tier_filter="MODERATE",
                                  out_path=MARKERS_DIR / "tier_II_markers.fna")
    # Tier III markers
    meta_iii = write_marker_fasta(regions, ref_seqs,
                                  tier_filter="DIVERGED",
                                  out_path=MARKERS_DIR / "tier_III_markers.fna")
    # Combined
    meta_all = write_marker_fasta(regions, ref_seqs,
                                  tier_filter=None,
                                  out_path=MARKERS_DIR / "all_discriminatory.fna")

    # Metadata TSV
    if meta_all:
        meta_path = MARKERS_DIR / "marker_metadata.tsv"
        with open(meta_path, "w", newline="") as fh:
            writer = csv.DictWriter(fh, fieldnames=meta_all[0].keys(), delimiter="\t")
            writer.writeheader()
            writer.writerows(meta_all)
        log.info("Marker metadata: %s (%d rows)", meta_path, len(meta_all))

    log.info("Tier II: %d markers | Tier III: %d markers | Total: %d",
             len(meta_ii), len(meta_iii), len(meta_all))


if __name__ == "__main__":
    main()
