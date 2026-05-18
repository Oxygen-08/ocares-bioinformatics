#!/usr/bin/env python3
"""
Triangulation experiment suite — run all pre-submission validation tests:
  Exp 1: Regularized XGBoost — does tighter regularization close the train/test gap?
  Exp 2: K-12 relabelling — does independent biological label change AUROC materially?
  Exp 3: Spike-in sensitivity — how does per-tier AUROC vary across low/mid/high abundance?
  Exp 4: Label shuffle null control — confirms current AUROC is above chance.

Outputs:
  data/results/ml/regularized_cv_results.tsv
  data/results/ml/k12_label_cv_results.tsv
  data/results/ml/experiment_summary.tsv
  data/results/figures/experiment_comparison.png
"""

import logging
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import xgboost as xgb
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from sklearn.model_selection import StratifiedKFold, cross_validate
from sklearn.metrics import roc_auc_score, average_precision_score
from sklearn.utils import shuffle as sk_shuffle

warnings.filterwarnings("ignore")
logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s")
log = logging.getLogger(__name__)

BASE    = Path(__file__).parent.parent.parent
RESULTS = BASE / "data" / "results" / "ml"
FIGS    = BASE / "data" / "results" / "figures"
BIOSANITY = BASE / "data" / "results" / "bio_sanity"
BLAST_DIR = BASE / "data" / "results" / "blast_screen"
META_DIR  = BASE / "data" / "metagenome"

FEATURE_COLS = [
    "blastn_identity", "cai_score", "gc_delta", "srna_density",
    "align_coverage", "kmer_deviation", "presence_pathogenic",
    "presence_non_pathogenic", "pangenome_score", "anvio_cluster_score",
]
N_FOLDS = 5


# ── Helpers ─────────────────────────────────────────────────────────────────

def cv_xgb(X, y, params: dict, label: str) -> dict:
    """Run StratifiedKFold CV and return summary stats."""
    n_pos, n_neg = y.sum(), (1 - y).sum()
    spw = n_neg / max(n_pos, 1)
    base = dict(
        n_estimators=200, learning_rate=0.05, eval_metric="aucpr",
        scale_pos_weight=spw, random_state=42, verbosity=0,
        use_label_encoder=False,
    )
    base.update(params)
    clf = xgb.XGBClassifier(**base)
    skf = StratifiedKFold(n_splits=N_FOLDS, shuffle=True, random_state=42)

    train_aucs, test_aucs, train_aprs, test_aprs, train_f1s, test_f1s = [], [], [], [], [], []

    for fold, (tr, te) in enumerate(skf.split(X, y)):
        Xtr, Xte = X[tr], X[te]
        ytr, yte = y[tr], y[te]
        clf.fit(Xtr, ytr)
        tr_prob = clf.predict_proba(Xtr)[:, 1]
        te_prob = clf.predict_proba(Xte)[:, 1]
        train_aucs.append(roc_auc_score(ytr, tr_prob))
        test_aucs.append(roc_auc_score(yte, te_prob))
        train_aprs.append(average_precision_score(ytr, tr_prob))
        test_aprs.append(average_precision_score(yte, te_prob))
        # F1 at 0.5 threshold
        yte_pred = (te_prob >= 0.5).astype(int)
        tp = ((yte_pred == 1) & (yte == 1)).sum()
        fp = ((yte_pred == 1) & (yte == 0)).sum()
        fn = ((yte_pred == 0) & (yte == 1)).sum()
        prec = tp / max(tp + fp, 1)
        rec  = tp / max(tp + fn, 1)
        f1   = 2 * prec * rec / max(prec + rec, 1e-9)
        test_f1s.append(f1)

    return {
        "label": label,
        "n_pos": int(n_pos), "n_neg": int(n_neg),
        "train_auroc_mean": np.mean(train_aucs),
        "train_auroc_std":  np.std(train_aucs),
        "test_auroc_mean":  np.mean(test_aucs),
        "test_auroc_std":   np.std(test_aucs),
        "train_auprc_mean": np.mean(train_aprs),
        "test_auprc_mean":  np.mean(test_aprs),
        "test_f1_mean":     np.mean(test_f1s),
        "overfit_gap":      np.mean(train_aucs) - np.mean(test_aucs),
    }


# ── Experiment 1: Regularized XGBoost ───────────────────────────────────────

def exp1_regularized(df: pd.DataFrame) -> dict:
    log.info("Exp 1: Regularized XGBoost (max_depth=3, reg_lambda=1, min_child_weight=5)")
    X = df[FEATURE_COLS].values.astype(float)
    y = df["label"].values.astype(int)
    res = cv_xgb(X, y, {"max_depth": 3, "reg_lambda": 1.0, "min_child_weight": 5},
                 "Regularized (depth=3, L2=1)")
    log.info(f"  Train AUROC: {res['train_auroc_mean']:.3f} ± {res['train_auroc_std']:.3f}")
    log.info(f"  Test  AUROC: {res['test_auroc_mean']:.3f} ± {res['test_auroc_std']:.3f}")
    log.info(f"  Overfit gap: {res['overfit_gap']:.3f}")
    return res


# ── Experiment 2: K-12 Relabelling ──────────────────────────────────────────

def exp2_k12_relabel(df: pd.DataFrame) -> tuple:
    log.info("Exp 2: K-12 absence relabelling (positive = DIVERGED AND absent from all K-12)")
    sanity_file = BIOSANITY / "diverged_marker_summary.tsv"
    if not sanity_file.exists():
        log.warning(f"  K-12 sanity file not found: {sanity_file}. Skipping.")
        return None, None

    k12 = pd.read_csv(sanity_file, sep="\t")
    log.info(f"  Loaded K-12 absence table: {len(k12)} DIVERGED markers")

    # Merge on marker_id
    merged = df.merge(k12[["marker_id", "absent_from_k12"]], on="marker_id", how="left")
    # New label: 1 only if DIVERGED AND absent from K-12; 0 otherwise
    merged["label_k12"] = (
        (merged["tier"] == "DIVERGED") & (merged["absent_from_k12"] == True)
    ).astype(int)

    n_new_pos = merged["label_k12"].sum()
    n_orig_pos = df["label"].sum()
    log.info(f"  Original positive class: {n_orig_pos} (all DIVERGED)")
    log.info(f"  K-12 label positive class: {n_new_pos} (DIVERGED + K-12 absent)")
    log.info(f"  {n_orig_pos - n_new_pos} DIVERGED markers relabelled to negative (K-12 present)")

    X = merged[FEATURE_COLS].values.astype(float)
    y_orig = merged["label"].values.astype(int)
    y_k12  = merged["label_k12"].values.astype(int)

    # Original params (as in paper)
    res_orig = cv_xgb(X, y_orig, {"max_depth": 4}, "Original (DIVERGED = positive)")
    res_k12  = cv_xgb(X, y_k12,  {"max_depth": 4}, "K-12 absence label")

    log.info(f"  Original AUROC: {res_orig['test_auroc_mean']:.3f} ± {res_orig['test_auroc_std']:.3f}")
    log.info(f"  K-12 label AUROC: {res_k12['test_auroc_mean']:.3f} ± {res_k12['test_auroc_std']:.3f}")
    delta = res_k12["test_auroc_mean"] - res_orig["test_auroc_mean"]
    log.info(f"  Delta (K12 - Original): {delta:+.3f}")
    if abs(delta) < 0.03:
        log.info("  → AUROC stable under relabelling: label circularity is NOT inflating results")
    elif delta < 0:
        log.info("  → AUROC drops under K-12 label: some circularity present in original label")
    else:
        log.info("  → AUROC improves under K-12 label: biologically cleaner positive class")

    return res_orig, res_k12


# ── Experiment 3: Spike-in Sensitivity ──────────────────────────────────────

def exp3_spike_sensitivity() -> pd.DataFrame:
    log.info("Exp 3: Spike-in sensitivity — BLAST AUROC per abundance condition")
    read_labels_file = META_DIR / "read_labels.tsv"
    blast_file       = BLAST_DIR / "blast_diverged.tsv"

    if not read_labels_file.exists() or not blast_file.exists():
        log.warning("  Missing read_labels.tsv or blast_diverged.tsv. Skipping.")
        return pd.DataFrame()

    # Load read labels
    rl = pd.read_csv(read_labels_file, sep="\t")
    log.info(f"  Read labels: {len(rl)} reads, conditions: {rl['condition'].unique().tolist()}")

    # Load BLAST hits — col 1 = subject (marker), col 0 = query (read)
    blast_cols = ["marker", "read_id", "pident", "length", "mismatch", "gapopen",
                  "qstart", "qend", "sstart", "send", "evalue", "bitscore"]
    blast = pd.read_csv(blast_file, sep="\t", names=blast_cols, comment="#")
    log.info(f"  BLAST hits: {len(blast)}")

    # Clean read_id — strip /1 /2 pair suffix for matching
    blast["read_id_clean"] = blast["read_id"].str.replace(r"/[12]$", "", regex=True)
    rl["read_id_clean"]    = rl["read_id"].str.replace(r"/[12]$", "", regex=True)

    # Mark reads that had any BLAST hit to DIVERGED markers
    hit_reads = set(blast["read_id_clean"].unique())
    rl["blast_hit"] = rl["read_id_clean"].isin(hit_reads).astype(int)

    rows = []
    for cond in ["low", "mid", "high"]:
        sub = rl[rl["condition"] == cond].copy()
        if len(sub) == 0:
            continue
        n_pathogen = sub["is_pathogen"].sum()
        n_total    = len(sub)
        pct_path   = 100 * n_pathogen / n_total

        if sub["is_pathogen"].nunique() < 2:
            auroc = float("nan")
        else:
            auroc = roc_auc_score(sub["is_pathogen"], sub["blast_hit"])

        tp = ((sub["blast_hit"] == 1) & (sub["is_pathogen"] == 1)).sum()
        fp = ((sub["blast_hit"] == 1) & (sub["is_pathogen"] == 0)).sum()
        fn = ((sub["blast_hit"] == 0) & (sub["is_pathogen"] == 1)).sum()
        tn = ((sub["blast_hit"] == 0) & (sub["is_pathogen"] == 0)).sum()
        sens = tp / max(tp + fn, 1)
        spec = tn / max(tn + fp, 1)

        rows.append({
            "condition": cond,
            "pct_pathogen": round(pct_path, 1),
            "n_reads": n_total,
            "n_pathogen_reads": int(n_pathogen),
            "AUROC": round(auroc, 4) if auroc == auroc else "N/A",
            "sensitivity": round(sens, 4),
            "specificity": round(spec, 4),
        })
        log.info(f"  [{cond}] O157={pct_path:.0f}%  AUROC={auroc:.3f}  "
                 f"sens={sens:.3f}  spec={spec:.3f}")

    return pd.DataFrame(rows)


# ── Experiment 4: Label Shuffle Null Control ─────────────────────────────────

def exp4_null_control(df: pd.DataFrame) -> dict:
    log.info("Exp 4: Label shuffle null control (expected AUROC ≈ 0.5)")
    X = df[FEATURE_COLS].values.astype(float)
    y = sk_shuffle(df["label"].values.astype(int), random_state=42)
    res = cv_xgb(X, y, {"max_depth": 4}, "Null control (shuffled labels)")
    log.info(f"  Null AUROC: {res['test_auroc_mean']:.3f} ± {res['test_auroc_std']:.3f}")
    return res


# ── Plotting ─────────────────────────────────────────────────────────────────

def plot_summary(results: list, spike_df: pd.DataFrame):
    fig = plt.figure(figsize=(14, 10))
    fig.patch.set_facecolor('#f8f9fa')

    # Panel A: Train vs Test AUROC comparison
    ax1 = fig.add_subplot(2, 2, 1)
    labels = [r["label"] for r in results if r]
    train_means = [r["train_auroc_mean"] for r in results if r]
    test_means  = [r["test_auroc_mean"]  for r in results if r]
    train_stds  = [r["train_auroc_std"]  for r in results if r]
    test_stds   = [r["test_auroc_std"]   for r in results if r]

    x = np.arange(len(labels))
    w = 0.35
    ax1.bar(x - w/2, train_means, w, yerr=train_stds, label="Train AUROC",
            color="#c0392b", alpha=0.8, capsize=4)
    ax1.bar(x + w/2, test_means, w, yerr=test_stds, label="Test AUROC",
            color="#2980b9", alpha=0.8, capsize=4)
    ax1.axhline(0.5, color='gray', linestyle='--', lw=1, label="Random baseline")
    ax1.set_xticks(x)
    ax1.set_xticklabels(labels, rotation=15, ha='right', fontsize=8)
    ax1.set_ylabel("AUROC")
    ax1.set_title("A: Train vs Test AUROC — All Experiments", fontweight='bold', fontsize=10)
    ax1.legend(fontsize=8)
    ax1.set_ylim(0.4, 1.05)
    ax1.set_facecolor('#f8f9fa')

    # Panel B: Overfitting gap
    ax2 = fig.add_subplot(2, 2, 2)
    gaps = [r["overfit_gap"] for r in results if r]
    colors = ['#c0392b' if g > 0.2 else '#e67e22' if g > 0.1 else '#27ae60'
              for g in gaps]
    bars = ax2.bar(range(len(labels)), gaps, color=colors, alpha=0.85)
    ax2.axhline(0.1, color='#e67e22', linestyle='--', lw=1, label="Moderate overfit threshold")
    ax2.axhline(0.2, color='#c0392b', linestyle='--', lw=1, label="Severe overfit threshold")
    ax2.set_xticks(range(len(labels)))
    ax2.set_xticklabels(labels, rotation=15, ha='right', fontsize=8)
    ax2.set_ylabel("Train AUROC − Test AUROC")
    ax2.set_title("B: Overfitting Gap per Configuration", fontweight='bold', fontsize=10)
    ax2.legend(fontsize=8)
    ax2.set_facecolor('#f8f9fa')

    # Panel C: Spike-in sensitivity
    ax3 = fig.add_subplot(2, 2, 3)
    if not spike_df.empty and "AUROC" in spike_df.columns:
        spike_df2 = spike_df[spike_df["AUROC"] != "N/A"].copy()
        spike_df2["AUROC"] = spike_df2["AUROC"].astype(float)
        ax3.plot(spike_df2["pct_pathogen"], spike_df2["AUROC"],
                 'o-', color='#8e44ad', lw=2, markersize=8, label="AUROC")
        ax3.plot(spike_df2["pct_pathogen"], spike_df2["sensitivity"],
                 's--', color='#27ae60', lw=2, markersize=8, label="Sensitivity")
        ax3.plot(spike_df2["pct_pathogen"], spike_df2["specificity"],
                 '^--', color='#2980b9', lw=2, markersize=8, label="Specificity")
        ax3.axhline(0.5, color='gray', linestyle=':', lw=1)
        ax3.set_xlabel("O157:H7 abundance in community (%)")
        ax3.set_ylabel("Performance metric")
        ax3.set_title("C: Spike-in Sensitivity — BLAST DIVERGED Screen", fontweight='bold', fontsize=10)
        ax3.legend(fontsize=8)
        ax3.set_ylim(0, 1.05)
    else:
        ax3.text(0.5, 0.5, "Spike-in data\nnot available", ha='center', va='center',
                 transform=ax3.transAxes, fontsize=12)
    ax3.set_facecolor('#f8f9fa')

    # Panel D: AUPRC comparison
    ax4 = fig.add_subplot(2, 2, 4)
    auprc_vals = [r["test_auprc_mean"] for r in results if r]
    random_auprc = 0.263  # 109/415
    ax4.bar(range(len(labels)), auprc_vals, color='#16a085', alpha=0.85)
    ax4.axhline(random_auprc, color='gray', linestyle='--', lw=1.5,
                label=f"Random baseline = {random_auprc:.3f}")
    ax4.set_xticks(range(len(labels)))
    ax4.set_xticklabels(labels, rotation=15, ha='right', fontsize=8)
    ax4.set_ylabel("AUPRC (test)")
    ax4.set_title("D: AUPRC — All Configurations", fontweight='bold', fontsize=10)
    ax4.legend(fontsize=8)
    ax4.set_facecolor('#f8f9fa')

    plt.suptitle(
        "Pre-Submission Validation: Regularization · Label Analysis · Spike-in Sensitivity",
        fontsize=11, fontweight='bold', y=1.01
    )
    plt.tight_layout()
    out = FIGS / "experiment_comparison.png"
    plt.savefig(out, dpi=140, bbox_inches='tight')
    plt.close()
    log.info(f"  Figure saved: {out}")


# ── Main ────────────────────────────────────────────────────────────────────

def main():
    log.info("Loading feature matrix...")
    df = pd.read_csv(RESULTS / "feature_matrix.tsv", sep="\t")
    log.info(f"  {len(df)} markers loaded | positives: {df['label'].sum()} | negatives: {(df['label']==0).sum()}")

    # Original baseline (for comparison)
    log.info("\n=== ORIGINAL BASELINE (max_depth=4, no extra regularization) ===")
    X = df[FEATURE_COLS].values.astype(float)
    y = df["label"].values.astype(int)
    res_orig = cv_xgb(X, y, {"max_depth": 4}, "Original (depth=4, paper params)")
    log.info(f"  Train: {res_orig['train_auroc_mean']:.3f} | Test: {res_orig['test_auroc_mean']:.3f} | Gap: {res_orig['overfit_gap']:.3f}")

    log.info("\n=== EXP 1: REGULARIZATION ===")
    res_reg = exp1_regularized(df)

    log.info("\n=== EXP 2: K-12 RELABELLING ===")
    res_orig2, res_k12 = exp2_k12_relabel(df)

    log.info("\n=== EXP 3: SPIKE-IN SENSITIVITY ===")
    spike_df = exp3_spike_sensitivity()

    log.info("\n=== EXP 4: NULL CONTROL ===")
    res_null = exp4_null_control(df)

    # Compile summary
    all_results = [r for r in [res_orig, res_reg, res_k12, res_null] if r]
    summary = pd.DataFrame(all_results)
    out_summary = RESULTS / "experiment_summary.tsv"
    summary.to_csv(out_summary, sep="\t", index=False)
    log.info(f"\nSummary saved: {out_summary}")

    if not spike_df.empty:
        spike_df.to_csv(RESULTS / "spike_sensitivity.tsv", sep="\t", index=False)

    # Print final table
    print("\n" + "=" * 80)
    print("EXPERIMENT SUMMARY")
    print("=" * 80)
    cols = ["label", "n_pos", "test_auroc_mean", "test_auroc_std",
            "train_auroc_mean", "overfit_gap", "test_auprc_mean", "test_f1_mean"]
    print(summary[cols].to_string(index=False, float_format=lambda x: f"{x:.3f}"))

    if not spike_df.empty:
        print("\nSPIKE-IN SENSITIVITY:")
        print(spike_df.to_string(index=False))

    # Interpret
    print("\n" + "=" * 80)
    print("INTERPRETATION")
    print("=" * 80)
    reg_gap  = res_reg["overfit_gap"]
    orig_gap = res_orig["overfit_gap"]
    print(f"Regularization effect on overfit gap: {orig_gap:.3f} → {reg_gap:.3f} "
          f"({'reduced' if reg_gap < orig_gap else 'unchanged'})")
    if res_k12:
        delta = res_k12["test_auroc_mean"] - res_orig["test_auroc_mean"]
        print(f"K-12 relabelling AUROC delta: {delta:+.3f} "
              f"({'< 0.03 — circularity NOT inflating result' if abs(delta) < 0.03 else '≥ 0.03 — circularity present'})")
    print(f"Null control AUROC: {res_null['test_auroc_mean']:.3f} "
          f"(should be ≈ 0.5 — confirms true signal)")

    plot_summary(all_results, spike_df)
    log.info("\nAll experiments complete.")


if __name__ == "__main__":
    main()
