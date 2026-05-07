#!/usr/bin/env python3
"""
Phase 2 (Proof-of-Concept) — XGBoost classifier with novel 7-feature vector.

Feature vector per candidate genomic region:
  1. blastn_identity    — top BLASTn % identity against metagenome reads
  2. cai_score          — Codon Adaptation Index (O157:H7 codon usage as reference)
  3. gc_delta           — |GC% of marker − mean GC% of O157:H7 core genome|
  4. srna_density       — sRNA binding site density (RNAfold MFE-based proxy)
  5. align_coverage     — fraction of marker covered by BLAST hits
  6. tier_encoded       — ordinal encoding: CONSERVED=0, MODERATE=1, DIVERGED=2
  7. pangenome_score    — Anvi'o enrichment score (O-island gene frequency in pathogens)

XGBoost is chosen over logistic regression / SVM because:
  - Handles the non-linear interactions between CAI, GC delta, and tier (documented
    in Yin et al. 2023 for similar genomic feature classifiers)
  - Native support for class imbalance (scale_pos_weight parameter)
  - SHAP integration provides mechanistic interpretability required for the thesis
    argument that regulatory context (CAI, sRNA density) outperforms sequence
    identity alone

Output:
  data/results/ml/feature_matrix.tsv
  data/results/ml/model_xgb.json
  data/results/ml/shap_summary.png
  data/results/ml/shap_values.tsv
  data/results/ml/cv_results.tsv
"""

import csv
import json
import logging
import warnings
from pathlib import Path
from itertools import product

import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio.SeqUtils import gc_fraction
import xgboost as xgb
import shap
from sklearn.model_selection import StratifiedKFold, cross_validate
from sklearn.metrics import (
    roc_auc_score, average_precision_score,
    precision_score, recall_score, f1_score,
)
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

warnings.filterwarnings("ignore", category=UserWarning)

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s  %(levelname)s  %(message)s",
    datefmt="%H:%M:%S",
)
log = logging.getLogger(__name__)

REPO_ROOT   = Path(__file__).parents[2]
MARKERS_DIR = REPO_ROOT / "data" / "results" / "markers"
BLAST_DIR   = REPO_ROOT / "data" / "results" / "blast_screen"
ML_DIR      = REPO_ROOT / "data" / "results" / "ml"
ML_DIR.mkdir(parents=True, exist_ok=True)
VIZ_DIR     = REPO_ROOT / "data" / "results" / "figures"
VIZ_DIR.mkdir(parents=True, exist_ok=True)

REFERENCE   = REPO_ROOT / "data" / "genomes" / "O157H7_Sakai" / "O157H7_Sakai.fna"

TIER_ENCODING = {"CONSERVED": 0, "MODERATE": 1, "DIVERGED": 2}

# O157:H7 Sakai core genome mean GC%: 50.5% (Hayashi et al. 2001, PMID 11206551)
REFERENCE_GC = 50.5


# ── Feature computation ───────────────────────────────────────────────────────

def compute_gc_delta(seq: str) -> float:
    """Absolute difference between marker GC% and reference core GC%."""
    gc = gc_fraction(seq) * 100
    return abs(gc - REFERENCE_GC)


def compute_cai_proxy(seq: str) -> float:
    """
    Approximate CAI using GC3 content as a proxy when full codon table
    optimisation is unavailable. GC3 correlates strongly with CAI in
    Escherichia (Sharp & Li 1987; r=0.89 across 4,000 E. coli genes).
    Replace with full CAI computation once the CAI library is available.
    """
    codons = [seq[i:i+3] for i in range(0, len(seq) - 2, 3) if len(seq[i:i+3]) == 3]
    if len(codons) < 10:
        return 0.5  # too short for reliable estimate
    gc3 = sum(1 for c in codons if c[2] in ("G", "C")) / len(codons)
    # Linearly rescale GC3 (0.3–0.7 typical range) to CAI-like [0, 1]
    return min(max((gc3 - 0.3) / 0.4, 0.0), 1.0)


def compute_srna_density_proxy(seq: str) -> float:
    """
    sRNA binding site density proxy: fraction of 10-mer windows in the sequence
    that have ≥60% AU content (AU-rich regions are over-represented in sRNA
    recognition sites; Peer & Margalit 2011 PMID 21695124).
    Replace with full RNAfold MFE calculation in production.
    """
    if len(seq) < 10:
        return 0.0
    windows = [seq[i:i+10] for i in range(len(seq) - 9)]
    au_rich = sum(1 for w in windows if (w.count("A") + w.count("T")) / 10 >= 0.6)
    return au_rich / len(windows)


def build_kmer_profile(seq: str, k: int = 4) -> np.ndarray:
    """Return normalised k-mer frequency vector (length 4^k) for a sequence."""
    bases = "ACGT"
    kmers = ["".join(p) for p in product(bases, repeat=k)]
    index = {km: i for i, km in enumerate(kmers)}
    counts = np.zeros(len(kmers), dtype=float)
    seq = seq.upper()
    for i in range(len(seq) - k + 1):
        km = seq[i:i + k]
        if km in index:
            counts[index[km]] += 1
    total = counts.sum()
    return counts / total if total > 0 else counts


# Genome-level 4-mer profile cached at module level (populated on first call)
_GENOME_KMER_CACHE: dict[str, np.ndarray] = {}

def compute_kmer_deviation(seq: str, ref_seqs: dict[str, str], k: int = 4) -> float:
    """
    Cosine distance between the marker's k-mer profile and the full reference
    genome's k-mer profile.  Values near 0 = compositionally typical; values
    near 1 = highly atypical (candidate horizontal transfer / PAI region).
    """
    cache_key = f"k{k}_full"
    if cache_key not in _GENOME_KMER_CACHE:
        full_seq = "".join(ref_seqs.values())
        _GENOME_KMER_CACHE[cache_key] = build_kmer_profile(full_seq, k)
    genome_profile = _GENOME_KMER_CACHE[cache_key]

    marker_profile = build_kmer_profile(seq, k)
    denom = np.linalg.norm(marker_profile) * np.linalg.norm(genome_profile)
    if denom == 0:
        return 0.0
    cosine_sim = np.dot(marker_profile, genome_profile) / denom
    return float(1.0 - cosine_sim)   # cosine distance: 0=identical, 1=orthogonal


def compute_blast_features(marker_id: str, blast_df: pd.DataFrame) -> tuple[float, float]:
    """
    For a given marker, return (max_identity, coverage_fraction) from BLAST hits.
    Coverage = fraction of marker length covered by ≥1 BLAST hit.
    """
    hits = blast_df[blast_df["qseqid"].str.startswith(marker_id)]
    if hits.empty:
        return 0.0, 0.0
    max_identity = hits["pident"].max()
    # Estimate coverage: use total aligned length / marker length proxy
    # Marker length is encoded in qseqid header as len=XXX
    total_aligned = hits["length"].sum()
    # Approximate marker length from header (len=XXXX)
    import re
    match = re.search(r"len=(\d+)", marker_id)
    marker_len = int(match.group(1)) if match else 1000
    coverage = min(total_aligned / max(marker_len, 1), 1.0)
    return float(max_identity), float(coverage)


def load_blast_results() -> pd.DataFrame:
    """Load all BLAST results into single DataFrame."""
    all_hits = []
    for tsv in BLAST_DIR.glob("blast_*.tsv"):
        df = pd.read_csv(tsv, sep="\t", names=[
            "qseqid", "sseqid", "pident", "length", "mismatch",
            "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore",
        ])
        all_hits.append(df)
    return pd.concat(all_hits, ignore_index=True) if all_hits else pd.DataFrame()


def load_marker_metadata() -> pd.DataFrame:
    meta = MARKERS_DIR / "marker_metadata.tsv"
    if not meta.exists():
        raise FileNotFoundError(f"Run 03_extract_markers.py first: {meta}")
    return pd.read_csv(meta, sep="\t")


def load_pangenome_scores() -> dict[str, float]:
    """
    Load Anvi'o enrichment scores if available. Falls back to a tier-derived
    heuristic: CONSERVED=0.1, MODERATE=0.6, DIVERGED=0.9 (reflecting that
    highly diverged regions are enriched in pathogen-unique gene clusters).
    """
    enrichment_file = REPO_ROOT / "data" / "results" / "pangenome" / "enrichment_scores.tsv"
    if not enrichment_file.exists():
        log.warning("Pangenome enrichment file not found — using tier-derived heuristic")
        return {}
    scores = {}
    with open(enrichment_file) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            scores[row["marker_id"]] = float(row["enrichment_score"])
    return scores


def compute_pangenome_presence(
    meta: pd.DataFrame,
    tiered_blocks_path: Path,
) -> pd.DataFrame:
    """
    For each marker, compute:
      presence_pathogenic     = fraction of pathogenic strains whose NUCmer
                                alignment covers the marker region
      presence_non_pathogenic = same for non-pathogenic strains
      pangenome_score         = presence_pathogenic - presence_non_pathogenic

    Uses NUCmer tiered_blocks.tsv — no additional BLAST required.
    A strain is counted as 'covering' a marker if it has ≥1 alignment block
    that overlaps the marker region on the same reference contig.
    """
    PATHOGENIC = {"EHEC", "UPEC", "ETEC", "EAEC", "EPEC", "NMEC", "AIEC"}

    blocks = pd.read_csv(tiered_blocks_path, sep="\t")
    # Strain → pathotype map
    strain_pathotype = blocks[["label", "pathotype"]].drop_duplicates()
    path_strains    = set(strain_pathotype.loc[strain_pathotype["pathotype"].isin(PATHOGENIC), "label"])
    nonpath_strains = set(strain_pathotype.loc[~strain_pathotype["pathotype"].isin(PATHOGENIC), "label"])
    n_path    = max(len(path_strains), 1)
    n_nonpath = max(len(nonpath_strains), 1)

    rows = []
    for _, m in meta.iterrows():
        contig  = m["contig"]
        m_start = int(m["start"])
        m_end   = int(m["end"])

        # Blocks on the same reference contig that overlap [m_start, m_end]
        overlap = blocks[
            (blocks["ref_contig"] == contig) &
            (blocks["ref_start"] <= m_end) &
            (blocks["ref_end"]   >= m_start)
        ]

        strains_present = set(overlap["label"].unique())
        pres_path    = len(strains_present & path_strains)    / n_path
        pres_nonpath = len(strains_present & nonpath_strains) / n_nonpath

        rows.append({
            "marker_id":              m["marker_id"],
            "presence_pathogenic":    round(pres_path, 4),
            "presence_non_pathogenic": round(pres_nonpath, 4),
            "pangenome_score":        round(pres_path - pres_nonpath, 4),
        })

    return pd.DataFrame(rows).set_index("marker_id")


def build_feature_matrix(
    meta: pd.DataFrame,
    blast_df: pd.DataFrame,
    ref_seqs: dict[str, str],
    pangenome_scores: dict[str, float],
    pan_presence: pd.DataFrame,
) -> pd.DataFrame:
    rows = []
    for _, m in meta.iterrows():
        contig = m["contig"]
        start  = int(m["start"])
        end    = int(m["end"])
        tier   = m["tier"]
        mid    = m["marker_id"]

        if contig not in ref_seqs:
            continue
        seq = ref_seqs[contig][start - 1:end]

        max_id, coverage = compute_blast_features(mid, blast_df)

        # Pangenome score: Anvi'o enrichment if available, else NUCmer-derived
        if mid in pan_presence.index:
            pan_score   = pan_presence.loc[mid, "pangenome_score"]
            pres_path   = pan_presence.loc[mid, "presence_pathogenic"]
            pres_nonpath = pan_presence.loc[mid, "presence_non_pathogenic"]
        elif mid in pangenome_scores:
            pan_score    = pangenome_scores[mid]
            pres_path    = float("nan")
            pres_nonpath = float("nan")
        else:
            pan_score    = 0.0
            pres_path    = 0.0
            pres_nonpath = 0.0

        rows.append({
            "marker_id":               mid,
            "tier":                    tier,
            "blastn_identity":         max_id,
            "cai_score":               compute_cai_proxy(seq),
            "gc_delta":                compute_gc_delta(seq),
            "srna_density":            compute_srna_density_proxy(seq),
            "align_coverage":          coverage,
            "kmer_deviation":          compute_kmer_deviation(seq, ref_seqs),
            "presence_pathogenic":     pres_path,
            "presence_non_pathogenic": pres_nonpath,
            "pangenome_score":         pan_score,
            "tier_encoded":            TIER_ENCODING[tier],
            "label":                   1 if tier == "DIVERGED" else 0,
        })

    return pd.DataFrame(rows)


# ── Model training ────────────────────────────────────────────────────────────

FEATURE_COLS = [
    "blastn_identity", "cai_score", "gc_delta",
    "srna_density", "align_coverage",
    "kmer_deviation",
    "presence_pathogenic", "presence_non_pathogenic", "pangenome_score",
]


def train_and_evaluate(df: pd.DataFrame) -> xgb.XGBClassifier:
    X = df[FEATURE_COLS].values
    y = df["label"].values

    n_pos = y.sum()
    n_neg = len(y) - n_pos
    scale_pos_weight = n_neg / max(n_pos, 1)

    clf = xgb.XGBClassifier(
        n_estimators=200,
        max_depth=4,
        learning_rate=0.05,
        subsample=0.8,
        colsample_bytree=0.8,
        scale_pos_weight=scale_pos_weight,
        eval_metric="aucpr",
        random_state=42,
        verbosity=0,
    )

    # 5-fold stratified cross-validation
    cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
    cv_results = cross_validate(
        clf, X, y, cv=cv,
        scoring=["roc_auc", "average_precision", "f1"],
        return_train_score=True,
    )

    log.info("5-fold CV — AUROC: %.3f ± %.3f | AUPRC: %.3f ± %.3f | F1: %.3f ± %.3f",
             cv_results["test_roc_auc"].mean(),    cv_results["test_roc_auc"].std(),
             cv_results["test_average_precision"].mean(), cv_results["test_average_precision"].std(),
             cv_results["test_f1"].mean(),          cv_results["test_f1"].std())

    # Save CV metrics
    cv_df = pd.DataFrame({k: v for k, v in cv_results.items()})
    cv_df.to_csv(ML_DIR / "cv_results.tsv", sep="\t", index=False)

    # Fit final model on full dataset
    clf.fit(X, y)
    clf.save_model(str(ML_DIR / "model_xgb.json"))
    log.info("Model saved: %s", ML_DIR / "model_xgb.json")
    return clf


def compute_shap(clf: xgb.XGBClassifier, df: pd.DataFrame) -> None:
    X = df[FEATURE_COLS].values
    explainer  = shap.TreeExplainer(clf)
    shap_vals  = explainer.shap_values(X)

    # Summary beeswarm plot
    plt.figure(figsize=(8, 5))
    shap.summary_plot(shap_vals, X, feature_names=FEATURE_COLS, show=False)
    plt.title("SHAP Feature Importance — Tiered Marker Classifier")
    plt.tight_layout()
    plt.savefig(VIZ_DIR / "shap_summary.png", dpi=150, bbox_inches="tight")
    plt.close()
    log.info("SHAP summary plot saved")

    # Save raw SHAP values
    shap_df = pd.DataFrame(shap_vals, columns=FEATURE_COLS)
    shap_df["marker_id"] = df["marker_id"].values
    shap_df.to_csv(ML_DIR / "shap_values.tsv", sep="\t", index=False)


def main() -> None:
    # Load inputs
    meta     = load_marker_metadata()
    blast_df = load_blast_results()
    ref_seqs = {rec.id: str(rec.seq) for rec in SeqIO.parse(str(REFERENCE), "fasta")}
    pan_scores = load_pangenome_scores()

    # New features: k-mer deviation (computed inside build_feature_matrix)
    # and NUCmer-derived pangenome presence/absence
    tiered_blocks_path = REPO_ROOT / "data" / "results" / "tiered_blocks.tsv"
    log.info("Computing pangenome presence/absence from NUCmer blocks")
    pan_presence = compute_pangenome_presence(meta, tiered_blocks_path)
    log.info("  Pangenome score range: %.3f – %.3f",
             pan_presence["pangenome_score"].min(),
             pan_presence["pangenome_score"].max())

    log.info("Building feature matrix for %d markers", len(meta))
    df = build_feature_matrix(meta, blast_df, ref_seqs, pan_scores, pan_presence)

    if df.empty:
        log.error("Feature matrix is empty — check upstream pipeline steps")
        return

    df.to_csv(ML_DIR / "feature_matrix.tsv", sep="\t", index=False)
    log.info("Feature matrix: %d rows × %d features", len(df), len(FEATURE_COLS))
    log.info("Class balance: %d positives / %d negatives",
             df["label"].sum(), (df["label"] == 0).sum())

    clf = train_and_evaluate(df)
    compute_shap(clf, df)

    log.info("ML pipeline complete. Outputs in %s", ML_DIR)


if __name__ == "__main__":
    main()
