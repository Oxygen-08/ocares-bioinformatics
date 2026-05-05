#!/usr/bin/env python3
"""
Screen tiered discriminatory markers against a simulated metagenome via BLASTn,
then compute a confusion matrix and ROC/PR curves comparing the tiered-marker
approach against naive (all-markers-equal) classification.

Pipeline:
  1. Build BLASTn database from simulated metagenome reads (data/metagenome/)
  2. BLASTn each tier's marker FASTA against the database
  3. Parse hits → assign read-level pathogen calls using identity + coverage thresholds
  4. Compare calls to ground-truth labels from InSilicoSeq abundance file
  5. Compute confusion matrix + ROC / PR curves
  6. Write results to data/results/blast_screen/

Identity thresholds used for classification (matching the tiered framework):
  Tier I  hits (>=95%):  high confidence pathogen read
  Tier II hits (85-95%): moderate confidence — flagged, not definitive
  Tier III (< 85%):      low confidence — excluded from pathogen calls

This stratification is the mechanistic novelty of the thesis: rather than
using a flat identity cutoff (as Kraken2/MetaPhlAn4 do), we propagate
alignment tier information into the classification decision.
"""

import subprocess
import sys
import csv
import json
import logging
from pathlib import Path
from collections import defaultdict

BIN_DIR = Path(sys.executable).parent

import numpy as np
import pandas as pd
from sklearn.metrics import (
    confusion_matrix, roc_curve, precision_recall_curve,
    auc, average_precision_score, roc_auc_score,
)
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s  %(levelname)s  %(message)s",
    datefmt="%H:%M:%S",
)
log = logging.getLogger(__name__)

REPO_ROOT    = Path(__file__).parents[2]
METAGENOME   = REPO_ROOT / "data" / "metagenome"
MARKERS_DIR  = REPO_ROOT / "data" / "results" / "markers"
BLAST_DIR    = REPO_ROOT / "data" / "results" / "blast_screen"
BLAST_DIR.mkdir(parents=True, exist_ok=True)
VIZ_DIR      = REPO_ROOT / "data" / "results" / "figures"
VIZ_DIR.mkdir(parents=True, exist_ok=True)

# BLASTn output format 6 fields
BLAST_COLS = [
    "qseqid", "sseqid", "pident", "length", "mismatch",
    "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore",
]

IDENTITY_CUTOFFS = {
    "CONSERVED": 95.0,
    "MODERATE":  85.0,
    "DIVERGED":   0.0,
    "COMBINED":  85.0,  # combined set screened at moderate threshold
}
MIN_ALIGN_LEN = 100  # bp
EVALUE_CUTOFF = 1e-10


def build_blast_db(reads_fastas: list[Path]) -> Path:
    """Build nucleotide BLASTn database from all community metagenome FASTAs."""
    db_prefix  = BLAST_DIR / "metagenome_db"
    combined   = BLAST_DIR / "metagenome_combined.fasta"
    if (db_prefix.parent / (db_prefix.name + ".nhr")).exists():
        log.info("SKIP  BLASTn DB already built: %s", db_prefix)
        return db_prefix

    if not combined.exists():
        log.info("Combining %d community FASTAs into one DB input", len(reads_fastas))
        with open(combined, "wb") as out:
            for fasta in reads_fastas:
                with open(fasta, "rb") as f:
                    out.write(f.read())

    log.info("Building BLASTn database from combined metagenome")
    subprocess.run(
        [
            str(BIN_DIR / "makeblastdb"),
            "-in", str(combined),
            "-dbtype", "nucl",
            "-out", str(db_prefix),
            "-title", "simulated_metagenome",
        ],
        check=True, capture_output=True,
    )
    return db_prefix


def run_blastn(query_fasta: Path, db_prefix: Path, out_file: Path, threads: int = 4) -> Path:
    """Run BLASTn; return path to tabular output."""
    if out_file.exists():
        log.info("SKIP  BLASTn output exists: %s", out_file.name)
        return out_file

    log.info("BLASTn: %s vs metagenome DB", query_fasta.name)
    subprocess.run(
        [
            str(BIN_DIR / "blastn"),
            "-query", str(query_fasta),
            "-db", str(db_prefix),
            "-out", str(out_file),
            "-outfmt", "6",
            "-perc_identity", "70",   # broad net; tier filtering applied post-hoc
            "-evalue", str(EVALUE_CUTOFF),
            "-num_threads", str(threads),
            "-word_size", "11",
            "-max_target_seqs", "500",
        ],
        check=True, capture_output=True,
    )
    return out_file


def parse_blast(blast_file: Path) -> pd.DataFrame:
    """Parse BLAST outfmt6 into DataFrame; return empty DF if no hits."""
    if not blast_file.exists() or blast_file.stat().st_size == 0:
        return pd.DataFrame(columns=BLAST_COLS)
    return pd.read_csv(blast_file, sep="\t", names=BLAST_COLS, comment="#")


def assign_read_labels(
    hits: pd.DataFrame,
    tier: str,
    ground_truth: dict[str, int],
) -> tuple[list[int], list[int], list[float]]:
    """
    Evaluate classification over ALL labelled reads in ground_truth.
    Reads with a qualifying BLAST hit (pident >= cutoff, length >= MIN_ALIGN_LEN)
    are called pathogen (pred=1); all other labelled reads are called
    non-pathogen (pred=0), producing correct TN and FN counts.

    Returns parallel lists: y_true, y_pred, y_score (for ROC curves).
    """
    cutoff = IDENTITY_CUTOFFS[tier]

    # Build read_id → max qualifying hit identity; only reads WITH a qualifying hit are here
    hit_scores: dict[str, float] = {}
    if not hits.empty:
        filtered = hits[(hits["pident"] >= cutoff) & (hits["length"] >= MIN_ALIGN_LEN)]
        if not filtered.empty:
            best = filtered.groupby("sseqid")["pident"].max()
            hit_scores = best.to_dict()

    # Evaluate every labelled read; reads absent from hit_scores are called negative
    y_true, y_pred, y_score = [], [], []
    for read_id, true_label in ground_truth.items():
        if read_id in hit_scores:
            pred  = 1
            score = hit_scores[read_id] / 100.0
        else:
            pred  = 0
            score = 0.0
        y_true.append(true_label)
        y_pred.append(pred)
        y_score.append(score)

    return y_true, y_pred, y_score


def build_contig_genome_map(genome_dir: Path, pathogen_labels: set[str]) -> dict[str, int]:
    """
    Map FASTA contig IDs → is_pathogen by reading each genome's FASTA headers.
    ISS encodes the contig ID as the prefix of every read header.
    """
    import re
    contig_map: dict[str, int] = {}
    manifest = genome_dir / "genome_manifest.tsv"
    if not manifest.exists():
        return contig_map
    with open(manifest) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            if row["fasta_path"] == "FAILED":
                continue
            label = row["label"]
            is_path = 1 if label in pathogen_labels else 0
            fasta = Path(row["fasta_path"])
            if not fasta.exists():
                continue
            with open(fasta) as fa:
                for line in fa:
                    if line.startswith(">"):
                        contig_id = line[1:].split()[0]
                        contig_map[contig_id] = is_path
    log.info("Contig→genome map: %d contigs (%d pathogen)",
             len(contig_map), sum(contig_map.values()))
    return contig_map


def load_ground_truth(read_labels_tsv: Path, contig_map: dict[str, int]) -> dict[str, int]:
    """
    Build read_id → is_pathogen using contig_map.
    ISS read headers: {contig_id}_{fragment_start}_{read_n}/{pair}
    Strip trailing _{int}_{int}/{int} to recover the contig ID.
    """
    import re
    _strip = re.compile(r'_\d+_\d+/\d+$')
    truth: dict[str, int] = {}
    if not read_labels_tsv.exists():
        log.warning("read_labels.tsv not found: %s", read_labels_tsv)
        return truth
    with open(read_labels_tsv) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            read_id  = row["read_id"]
            contig   = _strip.sub("", read_id)
            truth[read_id] = contig_map.get(contig, 0)
    log.info("Ground truth: %d labelled reads (%d pathogen)",
             len(truth), sum(truth.values()))
    return truth


def plot_roc(results: dict, out_path: Path) -> None:
    """Plot ROC curves for each tier side by side."""
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    for ax, metric in zip(axes, ["roc", "pr"]):
        for tier, data in results.items():
            if not data["y_true"]:
                continue
            y_true  = np.array(data["y_true"])
            y_score = np.array(data["y_score"])
            if metric == "roc":
                fpr, tpr, _ = roc_curve(y_true, y_score)
                auroc = roc_auc_score(y_true, y_score)
                ax.plot(fpr, tpr, label=f"{tier} (AUC={auroc:.3f})")
                ax.set_xlabel("False Positive Rate")
                ax.set_ylabel("True Positive Rate")
                ax.set_title("ROC Curves — Tiered Marker Classification")
            else:
                prec, rec, _ = precision_recall_curve(y_true, y_score)
                ap = average_precision_score(y_true, y_score)
                ax.plot(rec, prec, label=f"{tier} (AP={ap:.3f})")
                ax.set_xlabel("Recall")
                ax.set_ylabel("Precision")
                ax.set_title("Precision-Recall Curves")
            ax.legend(fontsize=8)

    plt.tight_layout()
    plt.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close()
    log.info("ROC/PR figure saved: %s", out_path.name)


def plot_confusion(cm: np.ndarray, tier: str, out_path: Path) -> None:
    plt.figure(figsize=(5, 4))
    sns.heatmap(
        cm, annot=True, fmt="d", cmap="Blues",
        xticklabels=["Non-pathogen", "Pathogen"],
        yticklabels=["Non-pathogen", "Pathogen"],
    )
    plt.title(f"Confusion Matrix — {tier} markers")
    plt.ylabel("True label")
    plt.xlabel("Predicted label")
    plt.tight_layout()
    plt.savefig(out_path, dpi=150)
    plt.close()


def main() -> None:
    # Locate reads FASTA — InSilicoSeq writes <prefix>_R1.fastq + R2
    # Convert to FASTA if needed, or use merged file
    reads_files = sorted(METAGENOME.glob("*.fasta")) + sorted(METAGENOME.glob("*.fna"))
    if not reads_files:
        log.error(
            "No metagenome FASTA found in %s — run 05_simulate_metagenome.py first",
            METAGENOME,
        )
        return

    db_prefix = build_blast_db(reads_files)

    # Ground truth — map contig IDs → genome labels → is_pathogen
    PATHOGEN_LABELS = {"O157H7_Sakai"}
    genome_dir  = REPO_ROOT / "data" / "genomes"
    contig_map  = build_contig_genome_map(genome_dir, PATHOGEN_LABELS)
    gt_file     = METAGENOME / "read_labels.tsv"
    gt          = load_ground_truth(gt_file, contig_map)

    tiers = {
        "MODERATE":  MARKERS_DIR / "tier_II_markers.fna",
        "DIVERGED":  MARKERS_DIR / "tier_III_markers.fna",
        "COMBINED":  MARKERS_DIR / "all_discriminatory.fna",
    }

    results = {}
    for tier, marker_fasta in tiers.items():
        if not marker_fasta.exists():
            log.warning("Marker FASTA missing: %s", marker_fasta.name)
            continue

        blast_out = BLAST_DIR / f"blast_{tier.lower()}.tsv"
        run_blastn(marker_fasta, db_prefix, blast_out)
        hits = parse_blast(blast_out)
        log.info("Tier %s: %d BLAST hits", tier, len(hits))

        y_true, y_pred, y_score = assign_read_labels(hits, tier, gt)
        results[tier] = {
            "y_true": y_true, "y_pred": y_pred, "y_score": y_score,
            "n_hits": len(hits),
        }

        if y_true:
            cm = confusion_matrix(y_true, y_pred)
            plot_confusion(cm, tier, VIZ_DIR / f"confusion_{tier.lower()}.png")
            log.info(
                "  %s — TN=%d FP=%d FN=%d TP=%d",
                tier, cm[0, 0], cm[0, 1], cm[1, 0], cm[1, 1],
            )

    if any(r["y_true"] for r in results.values()):
        plot_roc(results, VIZ_DIR / "roc_pr_curves.png")

    # Write numeric summary
    summary = []
    for tier, data in results.items():
        if not data["y_true"]:
            continue
        y_true  = np.array(data["y_true"])
        y_pred  = np.array(data["y_pred"])
        y_score = np.array(data["y_score"])
        cm      = confusion_matrix(y_true, y_pred)
        tn, fp, fn, tp = cm.ravel()
        summary.append({
            "tier":      tier,
            "n_hits":    data["n_hits"],
            "TP": int(tp), "FP": int(fp), "TN": int(tn), "FN": int(fn),
            "sensitivity": round(tp / max(tp + fn, 1), 4),
            "specificity": round(tn / max(tn + fp, 1), 4),
            "precision":   round(tp / max(tp + fp, 1), 4),
            "AUROC": round(roc_auc_score(y_true, y_score), 4) if len(set(y_true)) > 1 else "N/A",
        })

    if summary:
        out = BLAST_DIR / "classification_metrics.tsv"
        with open(out, "w", newline="") as fh:
            writer = csv.DictWriter(fh, fieldnames=summary[0].keys(), delimiter="\t")
            writer.writeheader()
            writer.writerows(summary)
        log.info("Classification metrics: %s", out)
        for row in summary:
            log.info(
                "  %-10s  Sens=%.3f  Spec=%.3f  Prec=%.3f  AUROC=%s",
                row["tier"], row["sensitivity"], row["specificity"],
                row["precision"], row["AUROC"],
            )


if __name__ == "__main__":
    main()
