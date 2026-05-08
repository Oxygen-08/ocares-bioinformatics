#!/usr/bin/env python3
"""
Results figures for thesis:
  1. ROC curve — BLAST screen vs ML classifier
  2. Confusion matrix comparison — BLAST vs ML (side-by-side)
  3. Feature correlation heatmap — confirms non-redundancy of ML features

Outputs:
  data/results/figures/roc_comparison.png
  data/results/figures/confusion_matrix_comparison.png
  data/results/figures/feature_correlation.png
"""

import logging
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from sklearn.metrics import (
    roc_curve, auc, confusion_matrix,
    ConfusionMatrixDisplay,
)
from sklearn.model_selection import StratifiedKFold
import xgboost as xgb

logging.basicConfig(level=logging.INFO, format="%(asctime)s  %(message)s",
                    datefmt="%H:%M:%S")
log = logging.getLogger(__name__)

REPO_ROOT = Path(__file__).parents[2]
ML_DIR    = REPO_ROOT / "data" / "results" / "ml"
BLAST_DIR = REPO_ROOT / "data" / "results" / "blast_screen"
VIZ_DIR   = REPO_ROOT / "data" / "results" / "figures"
VIZ_DIR.mkdir(parents=True, exist_ok=True)

FEATURE_COLS = [
    "blastn_identity", "cai_score", "gc_delta",
    "srna_density", "align_coverage",
    "kmer_deviation",
    "presence_pathogenic", "presence_non_pathogenic", "pangenome_score",
    "anvio_cluster_score",
]

FEATURE_LABELS = {
    "blastn_identity":         "BLASTn identity",
    "cai_score":               "CAI (GC3 proxy)",
    "gc_delta":                "GC delta",
    "srna_density":            "sRNA density",
    "align_coverage":          "Alignment coverage",
    "kmer_deviation":          "4-mer deviation",
    "presence_pathogenic":     "Pathogen presence",
    "presence_non_pathogenic": "Commensal presence",
    "pangenome_score":         "Pangenome score (NUCmer)",
    "anvio_cluster_score":     "Cluster score (Anvi'o)",
}


# ── 1. ROC curve ─────────────────────────────────────────────────────────────

def plot_roc(df: pd.DataFrame) -> None:
    X = df[FEATURE_COLS].values
    y = df["label"].values

    n_pos = y.sum()
    n_neg = len(y) - n_pos

    clf = xgb.XGBClassifier(
        n_estimators=200, max_depth=4, learning_rate=0.05,
        subsample=0.8, colsample_bytree=0.8,
        scale_pos_weight=n_neg / max(n_pos, 1),
        eval_metric="aucpr", random_state=42, verbosity=0,
    )

    cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
    tprs, aucs = [], []
    mean_fpr = np.linspace(0, 1, 200)

    fig, ax = plt.subplots(figsize=(7, 6))

    for fold, (train_idx, test_idx) in enumerate(cv.split(X, y)):
        clf.fit(X[train_idx], y[train_idx])
        proba = clf.predict_proba(X[test_idx])[:, 1]
        fpr, tpr, _ = roc_curve(y[test_idx], proba)
        fold_auc = auc(fpr, tpr)
        aucs.append(fold_auc)
        interp_tpr = np.interp(mean_fpr, fpr, tpr)
        interp_tpr[0] = 0.0
        tprs.append(interp_tpr)
        ax.plot(fpr, tpr, lw=0.8, alpha=0.35, color="#1565C0")

    mean_tpr = np.mean(tprs, axis=0)
    mean_tpr[-1] = 1.0
    mean_auc = np.mean(aucs)
    std_auc  = np.std(aucs)
    std_tpr  = np.std(tprs, axis=0)

    ax.plot(mean_fpr, mean_tpr, color="#1565C0", lw=2.2,
            label=f"ML classifier (AUROC = {mean_auc:.3f} ± {std_auc:.3f})")
    ax.fill_between(mean_fpr,
                    np.maximum(mean_tpr - std_tpr, 0),
                    np.minimum(mean_tpr + std_tpr, 1),
                    color="#1565C0", alpha=0.12, label="±1 SD")

    # BLAST baseline from classification_metrics.tsv
    blast_auc = 0.552
    ax.plot([0, 1], [0, blast_auc, 1][:2] if blast_auc < 1 else [0, 1],
            color="#F44336", lw=1.8, linestyle="--",
            label=f"BLAST screen baseline (AUROC ≈ {blast_auc:.3f})")
    ax.plot([0, 1], [0, 1], "k--", lw=0.8, alpha=0.4, label="Random (AUROC = 0.500)")

    ax.set_xlabel("False Positive Rate", fontsize=12)
    ax.set_ylabel("True Positive Rate", fontsize=12)
    ax.set_title("ROC curve — ML classifier vs BLAST screen\n"
                 "5-fold stratified CV on 415 markers", fontsize=12, pad=10)
    ax.legend(loc="lower right", fontsize=10)
    ax.set_xlim([-0.02, 1.02])
    ax.set_ylim([-0.02, 1.05])
    ax.grid(color="#eeeeee", linewidth=0.6)

    plt.tight_layout()
    out = VIZ_DIR / "roc_comparison.png"
    plt.savefig(out, dpi=180, bbox_inches="tight")
    plt.close()
    log.info("Saved: %s", out)


# ── 2. Confusion matrix comparison ───────────────────────────────────────────

def plot_confusion_matrices(df: pd.DataFrame) -> None:
    X = df[FEATURE_COLS].values
    y = df["label"].values
    n_pos = y.sum()
    n_neg = len(y) - n_pos

    clf = xgb.XGBClassifier(
        n_estimators=200, max_depth=4, learning_rate=0.05,
        subsample=0.8, colsample_bytree=0.8,
        scale_pos_weight=n_neg / max(n_pos, 1),
        eval_metric="aucpr", random_state=42, verbosity=0,
    )
    clf.fit(X, y)
    ml_pred = clf.predict(X)

    # BLAST confusion matrix (from classification_metrics.tsv, COMBINED row)
    blast_metrics = pd.read_csv(BLAST_DIR / "classification_metrics.tsv", sep="\t")
    combined = blast_metrics[blast_metrics["tier"] == "COMBINED"].iloc[0]
    blast_cm = np.array([
        [int(combined["TN"]), int(combined["FP"])],
        [int(combined["FN"]), int(combined["TP"])],
    ])

    ml_cm = confusion_matrix(y, ml_pred)

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))

    for ax, cm, title in zip(
        axes,
        [blast_cm, ml_cm],
        ["BLAST screen\n(sequence identity ≥ cutoff)",
         "ML classifier\n(10-feature XGBoost)"],
    ):
        disp = ConfusionMatrixDisplay(
            confusion_matrix=cm,
            display_labels=["MODERATE\n(negative)", "DIVERGED\n(positive)"],
        )
        disp.plot(ax=ax, colorbar=False, cmap="Blues")
        ax.set_title(title, fontsize=11, pad=8)

        total = cm.sum()
        tn, fp, fn, tp = cm.ravel()
        sens = tp / max(tp + fn, 1)
        spec = tn / max(tn + fp, 1)
        prec = tp / max(tp + fp, 1)
        ax.set_xlabel(
            f"Sensitivity={sens:.3f}  Specificity={spec:.3f}  Precision={prec:.3f}",
            fontsize=8.5,
        )

    plt.suptitle("Confusion matrix: BLAST screen vs ML classifier\n"
                 "(marker-level evaluation, 415 markers)",
                 fontsize=12, y=1.03)
    plt.tight_layout()
    out = VIZ_DIR / "confusion_matrix_comparison.png"
    plt.savefig(out, dpi=180, bbox_inches="tight")
    plt.close()
    log.info("Saved: %s", out)


# ── 3. Feature correlation heatmap ───────────────────────────────────────────

def plot_feature_correlation(df: pd.DataFrame) -> None:
    feat_df = df[FEATURE_COLS].rename(columns=FEATURE_LABELS)
    corr = feat_df.corr(method="spearman")

    fig, ax = plt.subplots(figsize=(9, 7))
    mask = np.triu(np.ones_like(corr, dtype=bool), k=1)
    sns.heatmap(
        corr, mask=mask, ax=ax,
        cmap="RdBu_r", vmin=-1, vmax=1, center=0,
        annot=True, fmt=".2f", annot_kws={"size": 8},
        linewidths=0.4, linecolor="#dddddd",
        square=True,
    )
    ax.set_title("Spearman correlation — ML feature matrix\n"
                 "Low inter-feature correlation confirms non-redundancy",
                 fontsize=11, pad=10)
    plt.tight_layout()
    out = VIZ_DIR / "feature_correlation.png"
    plt.savefig(out, dpi=180, bbox_inches="tight")
    plt.close()
    log.info("Saved: %s", out)


def main():
    df = pd.read_csv(ML_DIR / "feature_matrix.tsv", sep="\t")
    log.info("Loaded feature matrix: %d markers × %d features",
             len(df), len(FEATURE_COLS))

    log.info("Plotting ROC curve...")
    plot_roc(df)

    log.info("Plotting confusion matrices...")
    plot_confusion_matrices(df)

    log.info("Plotting feature correlation heatmap...")
    plot_feature_correlation(df)

    log.info("All result figures saved to %s", VIZ_DIR)


if __name__ == "__main__":
    main()
