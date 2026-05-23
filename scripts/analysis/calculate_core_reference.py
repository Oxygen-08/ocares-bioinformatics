#!/usr/bin/env python3
"""
Phase 1: Core Genome Profiling & Baseline Matrices
====================================================
Reads nucleotide sequences for all 1,962 CORE gene clusters from Anvi'o
contigs databases, then computes:

  1. CAI Reference Weight Matrix — relative adaptiveness w_ij per codon
     (Sharp & Li 1987: w_ij = f_ij / max_synonymous_freq_j)
  2. 59-dimensional Core RSCU Reference Vector
     (excluding stop codons, Met [ATG], Trp [TGG])

Output: data/metrics/core_evolutionary_reference.json

Run inside the fp_pipeline conda environment:
    conda run -n fp_pipeline python scripts/analysis/calculate_core_reference.py
"""

import json
import logging
import sqlite3
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
from Bio.Seq import Seq

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s  %(levelname)-7s  %(message)s",
    datefmt="%H:%M:%S",
)
log = logging.getLogger(__name__)

REPO_ROOT   = Path(__file__).parents[2]
PAN_DIR     = REPO_ROOT / "data" / "results" / "pangenome"
CONTIGS_DBS = PAN_DIR / "contigs_dbs"
GC_TSV      = PAN_DIR / "gene_clusters.tsv"
OUT_DIR     = REPO_ROOT / "data" / "metrics"
OUT_FILE    = OUT_DIR / "core_evolutionary_reference.json"

# ── Standard genetic code (DNA, sense strand) ────────────────────────────────

_GENETIC_CODE: dict[str, str] = {
    # Phe
    "TTT": "Phe", "TTC": "Phe",
    # Leu
    "TTA": "Leu", "TTG": "Leu",
    "CTT": "Leu", "CTC": "Leu", "CTA": "Leu", "CTG": "Leu",
    # Ile
    "ATT": "Ile", "ATC": "Ile", "ATA": "Ile",
    # Met — single synonymous codon, excluded from RSCU
    "ATG": "Met",
    # Val
    "GTT": "Val", "GTC": "Val", "GTA": "Val", "GTG": "Val",
    # Ser
    "TCT": "Ser", "TCC": "Ser", "TCA": "Ser", "TCG": "Ser",
    "AGT": "Ser", "AGC": "Ser",
    # Pro
    "CCT": "Pro", "CCC": "Pro", "CCA": "Pro", "CCG": "Pro",
    # Thr
    "ACT": "Thr", "ACC": "Thr", "ACA": "Thr", "ACG": "Thr",
    # Ala
    "GCT": "Ala", "GCC": "Ala", "GCA": "Ala", "GCG": "Ala",
    # Tyr
    "TAT": "Tyr", "TAC": "Tyr",
    # Stop codons — excluded
    "TAA": "Stop", "TAG": "Stop", "TGA": "Stop",
    # His
    "CAT": "His", "CAC": "His",
    # Gln
    "CAA": "Gln", "CAG": "Gln",
    # Asn
    "AAT": "Asn", "AAC": "Asn",
    # Lys
    "AAA": "Lys", "AAG": "Lys",
    # Asp
    "GAT": "Asp", "GAC": "Asp",
    # Glu
    "GAA": "Glu", "GAG": "Glu",
    # Cys
    "TGT": "Cys", "TGC": "Cys",
    # Trp — single synonymous codon, excluded from RSCU
    "TGG": "Trp",
    # Arg
    "CGT": "Arg", "CGC": "Arg", "CGA": "Arg", "CGG": "Arg",
    "AGA": "Arg", "AGG": "Arg",
    # Gly
    "GGT": "Gly", "GGC": "Gly", "GGA": "Gly", "GGG": "Gly",
}

# Codons excluded from the RSCU vector (stops + single-codon AAs)
_EXCLUDED: frozenset[str] = frozenset({"TAA", "TAG", "TGA", "ATG", "TGG"})

# Deterministic 59-codon RSCU vector order (alphabetical within standard genetic code)
RSCU_CODONS: list[str] = sorted(
    c for c, aa in _GENETIC_CODE.items()
    if c not in _EXCLUDED
)

# Synonymous groups for the 59 included codons
_SYN_GROUPS: dict[str, list[str]] = defaultdict(list)
for _c, _aa in _GENETIC_CODE.items():
    if _c not in _EXCLUDED:
        _SYN_GROUPS[_aa].append(_c)
_SYN_GROUPS = dict(_SYN_GROUPS)  # freeze


# ── Sequence utilities ────────────────────────────────────────────────────────

def reverse_complement(seq: str) -> str:
    return str(Seq(seq).reverse_complement())


def count_codons_from_cds(nt_seq: str, codon_counts: dict[str, int]) -> bool:
    """
    Tally sense-strand codons from a CDS (excluding the terminal stop codon).
    Returns True on success, False if sequence is malformed or too short.
    """
    if len(nt_seq) < 6 or len(nt_seq) % 3 != 0:
        return False
    for i in range(0, len(nt_seq) - 3, 3):  # -3 to skip stop codon
        codon = nt_seq[i:i + 3].upper()
        if len(codon) == 3 and codon.isalpha():
            codon_counts[codon] += 1
    return True


# ── Per-genome sequence extraction (batched SQL) ─────────────────────────────

def extract_genome_codon_counts(
    db_path: Path,
    gene_caller_ids: list[int],
    codon_counts: dict[str, int],
) -> tuple[int, int]:
    """
    Extract CDS sequences for a list of gene_caller_ids from an Anvi'o
    contigs database and accumulate codon counts.
    Returns (n_ok, n_skipped).
    """
    n_ok = n_skip = 0
    conn = sqlite3.connect(str(db_path))
    try:
        placeholders = ",".join("?" * len(gene_caller_ids))

        # Batch-fetch gene coordinates
        gene_rows = conn.execute(
            f"SELECT gene_callers_id, contig, start, stop, direction, partial "
            f"FROM genes_in_contigs WHERE gene_callers_id IN ({placeholders})",
            gene_caller_ids,
        ).fetchall()

        if not gene_rows:
            return 0, len(gene_caller_ids)

        # Collect required contigs
        needed_contigs = {row[1] for row in gene_rows if row[5] == 0}  # partial=0 only
        if not needed_contigs:
            return 0, len(gene_caller_ids)

        contig_placeholders = ",".join("?" * len(needed_contigs))
        contig_seqs = {
            row[0]: row[1]
            for row in conn.execute(
                f"SELECT contig, sequence FROM contig_sequences "
                f"WHERE contig IN ({contig_placeholders})",
                list(needed_contigs),
            ).fetchall()
        }

        for _, contig_name, start, stop, direction, partial in gene_rows:
            if partial:
                n_skip += 1
                continue
            if contig_name not in contig_seqs:
                n_skip += 1
                continue
            nt_seq = contig_seqs[contig_name][start:stop].upper()
            if direction == "r":
                nt_seq = reverse_complement(nt_seq)
            if count_codons_from_cds(nt_seq, codon_counts):
                n_ok += 1
            else:
                n_skip += 1
    finally:
        conn.close()
    return n_ok, n_skip


# ── RSCU and CAI computation ──────────────────────────────────────────────────

def compute_rscu(codon_counts: dict[str, int]) -> dict[str, float]:
    """
    RSCU_ij = observed_ij / (total_j / n_j)
    where total_j = sum of counts for all synonymous codons of amino acid j,
    and n_j = number of synonymous codons for amino acid j.

    Returns values only for the 59 RSCU_CODONS.
    """
    rscu: dict[str, float] = {}
    for aa, codons in _SYN_GROUPS.items():
        n = len(codons)
        total = sum(codon_counts.get(c, 0) for c in codons)
        if total == 0:
            for c in codons:
                rscu[c] = 1.0  # uniform fallback (no observed bias)
        else:
            expected_per_codon = total / n
            for c in codons:
                rscu[c] = codon_counts.get(c, 0) / expected_per_codon
    return rscu


def compute_cai_weights(codon_counts: dict[str, int]) -> dict[str, float]:
    """
    CAI relative adaptiveness (Sharp & Li 1987):
        w_ij = f_ij / max_k(f_kj)
    where f_ij is the frequency of codon i among synonymous codons for AA j.

    A pseudocount of 0.5 is applied per codon to prevent zero weights,
    which would cause the geometric mean to collapse to zero on rare codons.

    Returns values for the 59 non-excluded codons.
    """
    PSEUDOCOUNT = 0.5
    weights: dict[str, float] = {}
    for aa, codons in _SYN_GROUPS.items():
        freqs = {c: codon_counts.get(c, 0) + PSEUDOCOUNT for c in codons}
        max_freq = max(freqs.values())
        for c in codons:
            weights[c] = freqs[c] / max_freq
    return weights


# ── Main ─────────────────────────────────────────────────────────────────────

def main() -> None:
    log.info("=== Phase 1: Core Genome Profiling — fp_pipeline conda environment ===")
    log.info("RSCU vector dimensionality check: %d codons (expected 59)", len(RSCU_CODONS))
    assert len(RSCU_CODONS) == 59, f"RSCU vector is {len(RSCU_CODONS)}-dimensional; expected 59."

    log.info("Loading gene cluster assignments from %s", GC_TSV)
    gc = pd.read_csv(GC_TSV, sep="\t")

    n_total_genomes = gc["genome_name"].nunique()
    cluster_genome_counts = gc.groupby("gene_cluster_id")["genome_name"].nunique()
    core_ids = set(cluster_genome_counts[cluster_genome_counts == n_total_genomes].index)
    log.info(
        "Identified %d CORE clusters (present in all %d genomes) out of %d total",
        len(core_ids), n_total_genomes, gc["gene_cluster_id"].nunique(),
    )

    core_gc = gc[gc["gene_cluster_id"].isin(core_ids)].copy()
    log.info("Core gene-entry rows to process: %d", len(core_gc))

    codon_counts: dict[str, int] = defaultdict(int)
    total_ok = total_skip = 0
    genomes = sorted(core_gc["genome_name"].unique())

    for i, genome_name in enumerate(genomes, 1):
        # Anvi'o prefixes numeric genome names with "G_" in the TSV to avoid
        # column names starting with a digit; the DB files keep the original name.
        db_name = genome_name[2:] if genome_name.startswith("G_") and genome_name[2:].isdigit() else genome_name
        db_path = CONTIGS_DBS / f"{db_name}.db"
        if not db_path.exists():
            log.warning("[%d/%d] Missing contigs DB for %s — skipping genome",
                        i, len(genomes), genome_name)
            continue

        genome_rows = core_gc[core_gc["genome_name"] == genome_name]
        gene_ids = genome_rows["gene_caller_id"].astype(int).tolist()

        n_ok, n_skip = extract_genome_codon_counts(db_path, gene_ids, codon_counts)
        total_ok   += n_ok
        total_skip += n_skip

        log.info(
            "[%d/%d] %-25s  genes OK: %4d  skipped: %3d  running total: %d",
            i, len(genomes), genome_name, n_ok, n_skip, total_ok,
        )

    log.info("Extraction complete — sequences used: %d | skipped: %d", total_ok, total_skip)
    log.info("Total codons counted across core genome: %d", sum(codon_counts.values()))

    # ── Compute RSCU ──────────────────────────────────────────────────────────
    log.info("Computing 59-dimensional Core RSCU Reference Vector…")
    rscu_full = compute_rscu(codon_counts)
    rscu_vector = {c: round(rscu_full[c], 6) for c in RSCU_CODONS}

    rscu_vals = np.array(list(rscu_vector.values()))
    log.info(
        "RSCU stats — mean: %.4f  std: %.4f  min: %.4f  max: %.4f",
        rscu_vals.mean(), rscu_vals.std(), rscu_vals.min(), rscu_vals.max(),
    )

    # ── Compute CAI weights ───────────────────────────────────────────────────
    log.info("Computing CAI Reference Weight Matrix (Sharp & Li 1987)…")
    cai_weights = {c: round(v, 6) for c, v in compute_cai_weights(codon_counts).items()}

    w_vals = np.array(list(cai_weights.values()))
    log.info(
        "CAI weight stats — mean: %.4f  std: %.4f  min: %.4f  max: %.4f",
        w_vals.mean(), w_vals.std(), w_vals.min(), w_vals.max(),
    )

    # ── Save output ───────────────────────────────────────────────────────────
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    output = {
        "metadata": {
            "description": (
                "Core genome codon usage reference for E. coli O157:H7 61-genome pangenome. "
                "Used as baseline for Delta-CAI and RSCU Cosine Distance feature computation."
            ),
            "n_genomes": n_total_genomes,
            "n_core_clusters": len(core_ids),
            "n_sequences_used": total_ok,
            "n_sequences_skipped": total_skip,
            "total_codons_counted": sum(codon_counts.values()),
            "rscu_dimension": 59,
            "excluded_from_rscu": sorted(_EXCLUDED),
            "cai_pseudocount": 0.5,
            "reference": "Sharp & Li 1987 — doi:10.1093/nar/15.3.1281",
            "environment": "fp_pipeline conda env (Biopython, numpy, pandas)",
        },
        "rscu_codon_order": RSCU_CODONS,
        "rscu_reference_vector": rscu_vector,
        "cai_weight_matrix": cai_weights,
        "raw_codon_counts": {c: int(v) for c, v in codon_counts.items()},
    }

    with open(OUT_FILE, "w") as fh:
        json.dump(output, fh, indent=2)

    log.info("Reference profiles saved → %s", OUT_FILE)
    log.info("=== Phase 1 complete ===")


if __name__ == "__main__":
    main()
