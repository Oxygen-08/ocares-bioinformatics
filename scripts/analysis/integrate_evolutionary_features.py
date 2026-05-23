#!/usr/bin/env python3
"""
Phase 3: Matrix Integration — Cosine_Distance & Delta_CAI
==========================================================
Maps the per-cluster evolutionary features computed in Phase 2 onto the
existing marker-level ML feature matrix using the same Sakai gene overlap
logic used by compute_anvio_cluster_score() in 07_ml_classifier.py:

  For each NUCmer marker:
    1. Find all Sakai genes whose genomic coordinates overlap the marker.
    2. Look up each gene's Anvi'o gene cluster assignment.
    3. Retrieve the cluster's Cosine_Distance and Delta_CAI from Phase 2 output.
    4. Average across all overlapping clusters that have data.

Fill values for the two edge cases, documented here as the canonical record:
  • Core cluster genes overlapping a marker:
      Cosine_Distance = 0.0  (core is the reference by construction)
      Delta_CAI       = 0.0  (no deviation from core baseline)
  • Absent gene (no overlapping Sakai genes, or cluster not in Phase 2 output):
      Cosine_Distance = 1.0  (maximum penalty — represents maximum geometric divergence
                               from the core RSCU reference; gene is functionally absent
                               from this genomic context)
      Delta_CAI       = 1.0  (maximum penalty — represents total translational
                               incompatibility; absent gene contributes no expression
                               efficiency signal relative to the core genome baseline)

Modifies in-place:
  data/results/ml/feature_matrix.tsv
    → adds columns Cosine_Distance and Delta_CAI

Run inside the fp_pipeline conda environment:
    conda run -n fp_pipeline python scripts/analysis/integrate_evolutionary_features.py
"""

import logging
from pathlib import Path

import numpy as np
import pandas as pd
from Bio import SeqIO

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s  %(levelname)-7s  %(message)s",
    datefmt="%H:%M:%S",
)
log = logging.getLogger(__name__)

REPO_ROOT     = Path(__file__).parents[2]
PAN_DIR       = REPO_ROOT / "data" / "results" / "pangenome"
ML_DIR        = REPO_ROOT / "data" / "results" / "ml"
METRICS_DIR   = REPO_ROOT / "data" / "metrics"

GC_TSV          = PAN_DIR / "gene_clusters.tsv"
SAKAI_CALLS_TSV = PAN_DIR / "sakai_gene_calls.tsv"
REFERENCE       = REPO_ROOT / "data" / "genomes" / "O157H7_Sakai" / "O157H7_Sakai.fna"
REFERENCE_REFORM = PAN_DIR / "reformatted" / "O157H7_Sakai.fna"
FEATURES_TSV    = METRICS_DIR / "accessory_evolutionary_features.tsv"
MATRIX_TSV      = ML_DIR / "feature_matrix.tsv"

# Imputation constants for absent or unresolvable gene clusters.
# These values are applied when a marker has no overlapping Sakai genes
# or when a gene cluster has no Phase 2 feature data.
#
# Cosine_Distance = 1.0: maximum geometric penalty. A gene absent from the
#   genome contributes maximum codon-style divergence from the core reference —
#   the gene's translational machinery is entirely foreign or unrepresented.
#   Documented choice for peer review: absent ≡ maximally distant.
NEUTRAL_COS  = 1.0
#
# Delta_CAI = 1.0: maximum translational incompatibility penalty. An absent
#   gene contributes no CAI signal, and its absence is interpreted as total
#   deviation from core translation efficiency. This prevents the neutral
#   assumption (0.0) from masking the biological reality of gene absence.
#   Documented choice for peer review: absent ≡ maximally incompatible.
NEUTRAL_DCAI = 1.0
#
# Core cluster fill values: core genes ARE the reference by construction,
# so both features are identically zero (no distance, no deviation).
CORE_COS     = 0.0
CORE_DCAI    = 0.0


def build_contig_name_map() -> dict[str, str]:
    """
    Anvi'o reformats FASTA headers when building the contigs database.
    Reconstruct the reformatted → original mapping by FASTA order.
    """
    orig_ids     = [r.id for r in SeqIO.parse(str(REFERENCE),        "fasta")]
    reformed_ids = [r.id for r in SeqIO.parse(str(REFERENCE_REFORM), "fasta")]
    return dict(zip(reformed_ids, orig_ids))


def main() -> None:
    log.info("=== Phase 3: Matrix Integration — fp_pipeline conda environment ===")

    # ── Load Phase 2 output ───────────────────────────────────────────────────
    log.info("Loading per-cluster evolutionary features from %s", FEATURES_TSV)
    feat_df = pd.read_csv(FEATURES_TSV, sep="\t", index_col="gene_cluster_id")
    log.info("  %d accessory clusters with computed features", len(feat_df))

    # ── Load gene cluster assignments (Sakai genes only) ─────────────────────
    log.info("Loading gene cluster assignments…")
    gc = pd.read_csv(
        GC_TSV, sep="\t",
        usecols=["gene_caller_id", "gene_cluster_id", "genome_name"],
    )
    sakai_gc = (
        gc[gc["genome_name"] == "O157H7_Sakai"]
        [["gene_caller_id", "gene_cluster_id"]]
        .copy()
    )

    # Identify core cluster IDs (present in all genomes)
    n_total_genomes = gc["genome_name"].nunique()
    cluster_genome_counts = gc.groupby("gene_cluster_id")["genome_name"].nunique()
    core_cluster_ids = set(cluster_genome_counts[cluster_genome_counts == n_total_genomes].index)
    log.info("  Core clusters: %d | Accessory clusters in feature file: %d",
             len(core_cluster_ids), len(feat_df))

    # ── Load Sakai gene calls with genomic coordinates ────────────────────────
    log.info("Loading Sakai gene calls from %s", SAKAI_CALLS_TSV)
    gene_calls = pd.read_csv(
        SAKAI_CALLS_TSV, sep="\t",
        usecols=["gene_callers_id", "contig", "start", "stop"],
    )

    contig_map = build_contig_name_map()
    gene_calls["orig_contig"] = gene_calls["contig"].map(contig_map)
    gene_calls = gene_calls.dropna(subset=["orig_contig"])

    # Merge with cluster assignments
    gene_calls = gene_calls.merge(
        sakai_gc,
        left_on="gene_callers_id",
        right_on="gene_caller_id",
        how="left",
    )
    log.info("  Sakai gene calls merged: %d genes, %d with cluster assignment",
             len(gene_calls), gene_calls["gene_cluster_id"].notna().sum())

    # ── Load existing feature matrix ──────────────────────────────────────────
    log.info("Loading feature matrix from %s", MATRIX_TSV)
    matrix = pd.read_csv(MATRIX_TSV, sep="\t")
    log.info("  Feature matrix: %d markers × %d columns", *matrix.shape)

    if "Cosine_Distance" in matrix.columns or "Delta_CAI" in matrix.columns:
        log.info("  Columns Cosine_Distance / Delta_CAI already present — overwriting.")
        matrix = matrix.drop(columns=[c for c in ["Cosine_Distance", "Delta_CAI"]
                                       if c in matrix.columns])

    # ── Compute per-marker features ───────────────────────────────────────────
    cos_col:  list[float] = []
    dcai_col: list[float] = []

    n_core_only = n_missing = n_mixed = 0

    for _, row in matrix.iterrows():
        contig  = row["contig"]
        m_start = int(row["start"])
        m_end   = int(row["end"])

        # Sakai genes overlapping this marker
        overlap = gene_calls[
            (gene_calls["orig_contig"] == contig) &
            (gene_calls["start"]       <= m_end) &
            (gene_calls["stop"]        >= m_start)
        ]

        if overlap.empty:
            # IMPUTATION: no Sakai genes overlap this marker region.
            # Absent gene → maximum distance penalty for both features.
            # Cosine_Distance=1.0: maximum codon-style divergence from core.
            # Delta_CAI=1.0: maximum translational incompatibility from core.
            cos_col.append(NEUTRAL_COS)   # imputed: absent gene, max geometric divergence
            dcai_col.append(NEUTRAL_DCAI) # imputed: absent gene, max translational penalty
            n_missing += 1
            continue

        cos_vals:  list[float] = []
        dcai_vals: list[float] = []
        seen_core = False

        for _, gene in overlap.iterrows():
            cid = gene["gene_cluster_id"]
            if pd.isna(cid):
                continue
            if cid in core_cluster_ids:
                cos_vals.append(CORE_COS)
                dcai_vals.append(CORE_DCAI)
                seen_core = True
            elif cid in feat_df.index:
                cos_vals.append(feat_df.at[cid, "Cosine_Distance"])
                dcai_vals.append(feat_df.at[cid, "Delta_CAI"])
            # If cid not in feat_df and not core: skip (shouldn't occur normally)

        if not cos_vals:
            # IMPUTATION: Sakai genes overlapped but none had a resolvable
            # cluster ID (all NaN assignments). Treat as absent.
            # Cosine_Distance=1.0: maximum codon-style divergence from core.
            # Delta_CAI=1.0: maximum translational incompatibility from core.
            cos_col.append(NEUTRAL_COS)   # imputed: unresolvable cluster, max geometric divergence
            dcai_col.append(NEUTRAL_DCAI) # imputed: unresolvable cluster, max translational penalty
            n_missing += 1
        else:
            cos_col.append(float(np.mean(cos_vals)))
            dcai_col.append(float(np.mean(dcai_vals)))
            if seen_core and len(cos_vals) > 1:
                n_mixed += 1
            elif not seen_core:
                n_mixed += 1
            else:
                n_core_only += 1

    log.info(
        "Marker assignment summary: %d core-only | %d mixed/accessory | %d no-gene-data",
        n_core_only, n_mixed, n_missing,
    )

    # ── Append columns and save ───────────────────────────────────────────────
    matrix["Cosine_Distance"] = [round(v, 6) for v in cos_col]
    matrix["Delta_CAI"]       = [round(v, 6) for v in dcai_col]

    matrix.to_csv(MATRIX_TSV, sep="\t", index=False)
    log.info("Updated feature matrix saved → %s  (%d markers × %d columns)",
             MATRIX_TSV, *matrix.shape)

    # ── Phase 4: Statistical validation ───────────────────────────────────────
    log.info("=== Phase 4: Statistical Validation ===")

    feat_all = pd.read_csv(FEATURES_TSV, sep="\t")

    path_mask  = feat_all["pathogen_enriched"] == 1
    other_mask = ~path_mask

    def _report(series: pd.Series, label: str) -> None:
        log.info(
            "%-50s  n=%5d  mean=%+.4f  var=%.5f",
            label, len(series), series.mean(), series.var(),
        )

    _report(feat_all.loc[path_mask,  "Cosine_Distance"], "Cosine_Distance [pathogen-enriched]")
    _report(feat_all.loc[other_mask, "Cosine_Distance"], "Cosine_Distance [other accessory]")
    _report(feat_all.loc[path_mask,  "Delta_CAI"],       "Delta_CAI       [pathogen-enriched]")
    _report(feat_all.loc[other_mask, "Delta_CAI"],       "Delta_CAI       [other accessory]")

    log.info("=== Phase 3 + 4 complete ===")


if __name__ == "__main__":
    main()
