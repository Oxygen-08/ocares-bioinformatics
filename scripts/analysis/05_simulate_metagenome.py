#!/usr/bin/env python3
"""
Generate a simulated metagenomic community using InSilicoSeq 2.0.

Community design (justified against CAMI II strain-madness benchmark strategy):
  - O157:H7 Sakai          : target pathogen, spiked at 1%, 5%, 10% relative abundance
  - SE11 (EPEC)            : closest non-pathogenic relative (shares O-islands, no stx)
  - K-12 MG1655            : canonical non-pathogenic lab strain
  - Nissle 1917            : probiotic EPEC — deliberate false-positive trap
  - 2 gut commensals       : environmental background noise

Sequencing model: HiSeq 2500 (250 bp PE) — most common in CAMI II comparisons.
Three community compositions generated (low/mid/high spike) to characterise
sensitivity across a detection gradient.

Outputs:
  data/metagenome/community_low/     (O157:H7 at 1%)
  data/metagenome/community_mid/     (O157:H7 at 5%)
  data/metagenome/community_high/    (O157:H7 at 10%)
  data/metagenome/read_labels.tsv    (read_id → genome → is_pathogen)
"""

import subprocess
import csv
import logging
import random
import re
from pathlib import Path

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s  %(levelname)s  %(message)s",
    datefmt="%H:%M:%S",
)
log = logging.getLogger(__name__)

REPO_ROOT   = Path(__file__).parents[2]
GENOME_DIR  = REPO_ROOT / "data" / "genomes"
META_DIR    = REPO_ROOT / "data" / "metagenome"
META_DIR.mkdir(parents=True, exist_ok=True)

# Depth per community (total reads)
READ_DEPTHS = {
    "low":  500_000,
    "mid":  500_000,
    "high": 500_000,
}

# Community compositions — {strain_label: relative_abundance}
# Abundances sum to 1.0; O157:H7 varied; remainder distributed across background
COMMUNITIES = {
    "low": {
        "O157H7_Sakai":  0.01,
        "SE11":          0.30,
        "K12_MG1655":    0.30,
        "Nissle1917":    0.20,
        "HVM2062":       0.10,
        "HVM1157":       0.09,
    },
    "mid": {
        "O157H7_Sakai":  0.05,
        "SE11":          0.28,
        "K12_MG1655":    0.28,
        "Nissle1917":    0.19,
        "HVM2062":       0.10,
        "HVM1157":       0.10,
    },
    "high": {
        "O157H7_Sakai":  0.10,
        "SE11":          0.25,
        "K12_MG1655":    0.25,
        "Nissle1917":    0.20,
        "HVM2062":       0.10,
        "HVM1157":       0.10,
    },
}

PATHOGEN_STRAINS = {"O157H7_Sakai"}


def genome_path(label: str) -> Path:
    p = GENOME_DIR / label / f"{label}.fna"
    if not p.exists():
        raise FileNotFoundError(f"Genome not found: {p}")
    return p


def write_abundance_file(community: dict[str, float], out_file: Path) -> None:
    """Write InSilicoSeq abundance file: one line per genome — genome_path\tabundance."""
    with open(out_file, "w") as fh:
        for label, abund in community.items():
            fh.write(f"{genome_path(label)}\t{abund}\n")


def run_insilicoseq(
    abundance_file: Path,
    out_prefix: Path,
    n_reads: int,
    seed: int = 42,
) -> tuple[Path, Path]:
    """Run InSilicoSeq 2.0; return (R1_fastq, R2_fastq) paths."""
    r1 = out_prefix.parent / (out_prefix.name + "_R1.fastq")
    r2 = out_prefix.parent / (out_prefix.name + "_R2.fastq")

    if r1.exists() and r2.exists():
        log.info("SKIP  reads exist: %s", r1.parent.name)
        return r1, r2

    log.info("Simulating %d reads → %s", n_reads, out_prefix.parent.name)
    subprocess.run(
        [
            "iss", "generate",
            "--genomes", str(abundance_file),
            "--model", "HiSeq",
            "--output", str(out_prefix),
            "--n_reads", str(n_reads),
            "--cpus", "4",
            "--seed", str(seed),
        ],
        check=True,
    )
    return r1, r2


def fastq_to_fasta(fastq: Path, fasta: Path) -> Path:
    """Convert FASTQ to FASTA using awk (no Biopython dependency for large files)."""
    if fasta.exists():
        return fasta
    subprocess.run(
        f"awk 'NR%4==1{{print \">\"substr($0,2)}} NR%4==2{{print}}' {fastq} > {fasta}",
        shell=True, check=True,
    )
    return fasta


def build_read_label_table(
    communities: dict[str, dict[str, float]],
    meta_dir: Path,
) -> Path:
    """
    Parse InSilicoSeq read headers to extract genome origin, assign is_pathogen.
    InSilicoSeq 2.0 read headers: @<genome_file>_<read_number> (or similar).
    """
    out = meta_dir / "read_labels.tsv"
    if out.exists():
        return out

    rows = []
    for condition, community in communities.items():
        r1 = meta_dir / condition / f"{condition}_R1.fastq"
        if not r1.exists():
            continue
        with open(r1) as fh:
            for line in fh:
                if not line.startswith("@"):
                    continue
                read_id = line.strip().lstrip("@").split(" ")[0]
                # InSilicoSeq encodes source genome in header
                genome_label = "unknown"
                for label in community:
                    if label.lower().replace("_", "") in read_id.lower().replace("_", ""):
                        genome_label = label
                        break
                rows.append({
                    "read_id":    read_id,
                    "condition":  condition,
                    "genome":     genome_label,
                    "is_pathogen": 1 if genome_label in PATHOGEN_STRAINS else 0,
                })

    with open(out, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=["read_id", "condition", "genome", "is_pathogen"],
                                delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    log.info("Read label table: %s (%d reads labelled)", out, len(rows))
    return out


def main() -> None:
    # Verify genomes exist before starting
    for condition, community in COMMUNITIES.items():
        for label in community:
            try:
                genome_path(label)
            except FileNotFoundError as e:
                log.error("%s — run 01_download_genomes.py first", e)
                return

    for condition, community in COMMUNITIES.items():
        out_dir = META_DIR / condition
        out_dir.mkdir(exist_ok=True)

        abund_file = out_dir / "abundance.tsv"
        write_abundance_file(community, abund_file)

        n_reads = READ_DEPTHS[condition]
        prefix  = out_dir / condition

        try:
            r1, r2 = run_insilicoseq(abund_file, prefix, n_reads, seed=42)
        except (subprocess.CalledProcessError, FileNotFoundError) as exc:
            log.error("InSilicoSeq failed for %s: %s", condition, exc)
            continue

        # Merge R1+R2 into single FASTA for BLASTn
        merged_fastq = out_dir / f"{condition}_merged.fastq"
        merged_fasta = META_DIR / f"{condition}.fasta"
        if not merged_fasta.exists():
            subprocess.run(f"cat {r1} {r2} > {merged_fastq}", shell=True, check=True)
            fastq_to_fasta(merged_fastq, merged_fasta)
        log.info("Merged FASTA: %s", merged_fasta.name)

    build_read_label_table(COMMUNITIES, META_DIR)
    log.info("Simulation complete. Communities: %s", list(COMMUNITIES.keys()))


if __name__ == "__main__":
    main()
