"""
Kraken2 custom-database comparison for thesis Section 4.3.

Builds a minimal Kraken2 database from the 415 markers:
  - DIVERGED markers (n=109)  → taxid 2 (pathogen-specific)
  - MODERATE/CONSERVED (n=306) → taxid 3 (non-specific)

Classifies reads from each abundance condition (low=1%, mid=5%, high=10%).
Reports AUROC at read level vs is_pathogen ground truth from read_labels.tsv.
Compares to BLAST baseline (0.552) and XGBoost test AUROC (0.673).
"""

import os
import subprocess
import shutil
from pathlib import Path
import pandas as pd
import numpy as np
from sklearn.metrics import roc_auc_score, roc_curve, precision_recall_curve, auc
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

PROJ = Path("/Users/user/Documents/Projects/bioinformatics")
MARKERS_DIR  = PROJ / "data/results/markers"
META_DIR     = PROJ / "data/metagenome"
OUT_DIR      = PROJ / "data/results/kraken2"
FIG_DIR      = PROJ / "data/results/figures"
LABELS_TSV   = META_DIR / "read_labels.tsv"

ALL_MARKERS  = MARKERS_DIR / "all_discriminatory.fna"
TIER3        = MARKERS_DIR / "tier_III_markers.fna"    # 109 DIVERGED

DB_DIR       = OUT_DIR / "db"
KRAKEN2       = "kraken2"        # classification
KRAKEN2_BUILD = "kraken2-build"  # database construction

CONDITIONS = {
    "low":  META_DIR / "low/low_merged.fastq",
    "mid":  META_DIR / "mid/mid_merged.fastq",
    "high": META_DIR / "high/high_merged.fastq",
}

TAXID_PATHOGEN     = 2   # DIVERGED markers → O157-specific
TAXID_NON_PATHOGEN = 3   # MODERATE/CONSERVED → shared/core


# ── 1. Build minimal taxonomy ──────────────────────────────────────────────────

def build_taxonomy(db: Path) -> None:
    tax = db / "taxonomy"
    tax.mkdir(parents=True, exist_ok=True)

    # nodes.dmp: taxid | parent_taxid | rank | ...
    # Root(1) → pathogen(2) + nonpathogen(3)
    nodes = (
        "1\t|\t1\t|\tno rank\t|\t\t|\t8\t|\t0\t|\t1\t|\t0\t|\t0\t|\t0\t|\t0\t|\t0\t|\n"
        "2\t|\t1\t|\tspecies\t|\t\t|\t0\t|\t1\t|\t11\t|\t1\t|\t0\t|\t1\t|\t0\t|\t0\t|\n"
        "3\t|\t1\t|\tspecies\t|\t\t|\t0\t|\t1\t|\t11\t|\t1\t|\t0\t|\t1\t|\t0\t|\t0\t|\n"
    )
    names = (
        "1\t|\troot\t|\t\t|\tscientific name\t|\n"
        "2\t|\tEscherichia coli O157:H7\t|\t\t|\tscientific name\t|\n"
        "3\t|\tEscherichia coli commensal\t|\t\t|\tscientific name\t|\n"
    )
    (tax / "nodes.dmp").write_text(nodes)
    (tax / "names.dmp").write_text(names)
    # Empty accession2taxid (we tag via FASTA header)
    (tax / "seqid2taxid.map").write_text("")
    print("  Taxonomy written.")


# ── 2. Prepare library FASTA with kraken taxid headers ────────────────────────

def build_library_fasta(db: Path) -> Path:
    """
    Read all_discriminatory.fna; DIVERGED markers get |kraken:taxid|2 header,
    others get |kraken:taxid|3. Write to db/library/sequences.fna.
    """
    lib = db / "library"
    lib.mkdir(parents=True, exist_ok=True)
    out_fasta = lib / "sequences.fna"

    # Get DIVERGED marker IDs from tier_III_markers.fna headers
    diverged_ids = set()
    with open(TIER3) as f:
        for line in f:
            if line.startswith(">"):
                mid = line[1:].split("|")[0].strip()
                diverged_ids.add(mid)
    print(f"  DIVERGED marker IDs: {len(diverged_ids)}")

    n_pathogen = n_non_pathogen = 0
    with open(ALL_MARKERS) as fin, open(out_fasta, "w") as fout:
        current_taxid = None
        for line in fin:
            if line.startswith(">"):
                mid = line[1:].split("|")[0].strip()
                taxid = TAXID_PATHOGEN if mid in diverged_ids else TAXID_NON_PATHOGEN
                # Kraken2 requires: >seq_name|kraken:taxid|TAXID
                fout.write(f">{mid}|kraken:taxid|{taxid}\n")
                if taxid == TAXID_PATHOGEN:
                    n_pathogen += 1
                else:
                    n_non_pathogen += 1
            else:
                fout.write(line)
    print(f"  Library: {n_pathogen} pathogen seqs, {n_non_pathogen} non-pathogen seqs")
    return out_fasta


# ── 3. Build Kraken2 database ─────────────────────────────────────────────────

def kraken2_build(db: Path, lib_fasta: Path) -> None:
    # Step 1: add sequences to library
    cmd_add = [
        KRAKEN2_BUILD, "--db", str(db),
        "--add-to-library", str(lib_fasta),
    ]
    print(f"  Adding library: {' '.join(cmd_add)}")
    r = subprocess.run(cmd_add, capture_output=True, text=True)
    if r.returncode != 0:
        print(f"  STDERR:\n{r.stderr[-1000:]}")
        raise RuntimeError("kraken2-build --add-to-library failed")

    # Step 2: build the database
    cmd_build = [
        KRAKEN2_BUILD, "--db", str(db),
        "--build",
        "--kmer-len", "31",
        "--minimizer-len", "29",
        "--minimizer-spaces", "6",
        "--threads", "4",
    ]
    print(f"  Building DB: {' '.join(cmd_build)}")
    r = subprocess.run(cmd_build, capture_output=True, text=True)
    if r.returncode != 0:
        print(f"  STDERR:\n{r.stderr[-2000:]}")
        raise RuntimeError("kraken2-build --build failed")
    print("  DB built.")


# ── 4. Run Kraken2 on each condition ──────────────────────────────────────────

def run_kraken2(db: Path, fastq: Path, out_prefix: Path) -> Path:
    """Run Kraken2 on a merged FASTQ; return path to output classification file."""
    out_file   = Path(str(out_prefix) + "_classified.txt")
    report_file = Path(str(out_prefix) + "_report.txt")
    cmd = [
        KRAKEN2, "--db", str(db),
        "--output", str(out_file),
        "--report", str(report_file),
        "--threads", "4",
        str(fastq),
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"  STDERR:\n{result.stderr[-1000:]}")
        raise RuntimeError(f"kraken2 failed for {fastq.name}")
    # Print summary from stderr
    for line in result.stderr.strip().split("\n")[-5:]:
        print(f"    {line}")
    return out_file


# ── 5. Parse Kraken2 output → score per read ─────────────────────────────────

def parse_kraken2_output(classified_file: Path) -> pd.DataFrame:
    """
    Kraken2 output columns: C/U  read_id  taxid  length  kmer_string
    Score = fraction of k-mers assigned to pathogen taxid.
    """
    rows = []
    with open(classified_file) as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) < 5:
                continue
            status  = parts[0]   # C or U
            read_id = parts[1]
            taxid   = int(parts[2]) if parts[2].isdigit() else 0
            kmer_str = parts[4] if len(parts) > 4 else ""

            # Parse kmer hit breakdown: "2:15 3:8 0:5" → {taxid: count}
            kmer_counts = {}
            for chunk in kmer_str.split():
                if ":" in chunk:
                    t, c = chunk.split(":", 1)
                    try:
                        kmer_counts[int(t)] = kmer_counts.get(int(t), 0) + int(c)
                    except ValueError:
                        pass

            total_kmers = sum(kmer_counts.values())
            pathogen_kmers = kmer_counts.get(TAXID_PATHOGEN, 0)
            # Score = fraction of k-mers mapping to pathogen-specific taxid
            score = pathogen_kmers / total_kmers if total_kmers > 0 else 0.0
            rows.append({"read_id": read_id, "kraken_taxid": taxid,
                         "kraken_score": score, "is_classified": int(status == "C")})
    return pd.DataFrame(rows)


# ── 6. Evaluate per condition ─────────────────────────────────────────────────

def evaluate_condition(condition: str, db: Path) -> dict:
    fastq = CONDITIONS[condition]
    prefix = OUT_DIR / condition / condition
    (OUT_DIR / condition).mkdir(parents=True, exist_ok=True)

    print(f"\n  [{condition}] Running Kraken2 on {fastq.name}...")
    classified_file = run_kraken2(db, fastq, prefix)

    print(f"  [{condition}] Parsing output...")
    kraken_df = parse_kraken2_output(classified_file)

    # Load ground truth
    labels_df = pd.read_csv(LABELS_TSV, sep="\t")
    labels_cond = labels_df[labels_df["condition"] == condition][["read_id", "is_pathogen"]]

    # Merge — some reads may be absent from kraken output if completely unclassified
    merged = labels_cond.merge(kraken_df, on="read_id", how="left")
    merged["kraken_score"] = merged["kraken_score"].fillna(0.0)
    merged["is_classified"] = merged["is_classified"].fillna(0).astype(int)

    y_true = merged["is_pathogen"].values
    y_score = merged["kraken_score"].values

    # Binary classification AUROC
    try:
        auroc = roc_auc_score(y_true, y_score)
    except Exception:
        auroc = float("nan")

    # Sensitivity/specificity at binary threshold (classified by pathogen taxid)
    y_pred_binary = (merged["kraken_taxid"] == TAXID_PATHOGEN).astype(int)
    tp = ((y_pred_binary == 1) & (y_true == 1)).sum()
    fp = ((y_pred_binary == 1) & (y_true == 0)).sum()
    fn = ((y_pred_binary == 0) & (y_true == 1)).sum()
    tn = ((y_pred_binary == 0) & (y_true == 0)).sum()
    sensitivity = tp / (tp + fn) if (tp + fn) > 0 else 0.0
    specificity = tn / (tn + fp) if (tn + fp) > 0 else 0.0
    pct_pathogen = 100 * y_true.mean()

    result = {
        "condition": condition,
        "pct_pathogen": round(pct_pathogen, 1),
        "n_reads": len(merged),
        "n_pathogen_reads": int(y_true.sum()),
        "n_classified_pathogen": int(y_pred_binary.sum()),
        "AUROC": round(auroc, 4),
        "sensitivity": round(sensitivity, 4),
        "specificity": round(specificity, 4),
    }
    print(f"  [{condition}] O157={pct_pathogen:.0f}%  AUROC={auroc:.4f}  "
          f"sens={sensitivity:.3f}  spec={specificity:.3f}")
    return result


# ── 7. Plot comparison figure ─────────────────────────────────────────────────

def plot_comparison(results: list[dict]) -> None:
    fig, axes = plt.subplots(1, 3, figsize=(15, 5), facecolor="#0d1117")
    fig.suptitle("Kraken2 Custom-DB vs BLAST Baseline — Read-Level Detection",
                 color="white", fontsize=14, fontweight="bold")

    baselines = {
        "BLAST AUC": 0.552,
        "XGBoost AUC\n(marker-level)": 0.673,
    }
    conditions_order = ["low", "mid", "high"]

    for ax, result in zip(axes, [r for c in conditions_order for r in results if r["condition"] == c]):
        cond = result["condition"]
        pct = result["pct_pathogen"]
        auroc = result["AUROC"]
        sens = result["sensitivity"]
        spec = result["specificity"]

        bars = ax.bar(
            ["Kraken2\nAUROC", "BLAST\nAUROC", "XGBoost\nAUROC\n(marker)"],
            [auroc, 0.552, 0.673],
            color=["#58a6ff", "#f78166", "#3fb950"],
            alpha=0.85, width=0.5,
        )
        ax.axhline(0.5, color="white", linestyle="--", linewidth=0.8, alpha=0.5)
        ax.set_ylim(0, 1.0)
        ax.set_title(f"O157 = {pct:.0f}% spike-in", color="white", fontsize=11)
        ax.set_ylabel("AUROC", color="white")
        ax.tick_params(colors="white")
        ax.set_facecolor("#161b22")
        for spine in ax.spines.values():
            spine.set_edgecolor("#30363d")
        for bar, val in zip(bars, [auroc, 0.552, 0.673]):
            ax.text(bar.get_x() + bar.get_width() / 2, val + 0.02,
                    f"{val:.3f}", ha="center", va="bottom", color="white", fontsize=9)
        ax.text(0.5, -0.18, f"sens={sens:.2f}  spec={spec:.2f}",
                ha="center", transform=ax.transAxes, color="#8b949e", fontsize=8)

    plt.tight_layout()
    out_path = FIG_DIR / "kraken2_comparison.png"
    plt.savefig(out_path, dpi=150, bbox_inches="tight", facecolor=fig.get_facecolor())
    plt.close()
    print(f"\n  Figure saved: {out_path}")


# ── Main ──────────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    FIG_DIR.mkdir(parents=True, exist_ok=True)

    print("=== STEP 1: Build taxonomy ===")
    build_taxonomy(DB_DIR)

    print("\n=== STEP 2: Build library FASTA ===")
    build_library_fasta(DB_DIR)

    print("\n=== STEP 3: Build Kraken2 database ===")
    lib_fasta = DB_DIR / "library" / "sequences.fna"
    kraken2_build(DB_DIR, lib_fasta)

    print("\n=== STEP 4: Classify reads per condition ===")
    results = []
    for condition in ["low", "mid", "high"]:
        results.append(evaluate_condition(condition, DB_DIR))

    print("\n" + "=" * 70)
    print("KRAKEN2 COMPARISON RESULTS")
    print("=" * 70)
    df = pd.DataFrame(results)
    print(df.to_string(index=False))

    df.to_csv(OUT_DIR / "kraken2_results.tsv", sep="\t", index=False)
    print(f"\nResults saved: {OUT_DIR / 'kraken2_results.tsv'}")

    print("\n=== STEP 5: Plot ===")
    plot_comparison(results)

    print("\n" + "=" * 70)
    print("THESIS CONTEXT")
    print("=" * 70)
    mean_auroc = df["AUROC"].mean()
    print(f"Kraken2 mean AUROC (custom DB, 3 conditions): {mean_auroc:.3f}")
    print(f"BLAST baseline AUROC (from thesis):           0.552")
    print(f"XGBoost AUROC (marker-level, from thesis):    0.673")
    print(f"\nNote: Kraken2 and BLAST operate at read level (direct detection).")
    print(f"XGBoost operates at marker level (feature-based classification).")
    print(f"The three metrics measure different tasks but share the same ground truth.")
