#!/usr/bin/env python3
"""
Pangenome summary figures for Report 2:
  1. Gene cluster distribution — core / soft-core / accessory / unique
  2. Functional enrichment dot plot — PATHOGEN-enriched COG14 functions (q<0.05)
  3. Presence/absence heatmap — top 20 differentially present gene clusters

Outputs:
  data/results/figures/pangenome_composition.png
  data/results/figures/pangenome_enrichment.png
  data/results/figures/pangenome_heatmap.png
"""

import logging
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd
import seaborn as sns

logging.basicConfig(level=logging.INFO, format="%(asctime)s  %(message)s",
                    datefmt="%H:%M:%S")
log = logging.getLogger(__name__)

REPO_ROOT   = Path(__file__).parents[2]
PAN_DIR     = REPO_ROOT / "data" / "results" / "pangenome"
VIZ_DIR     = REPO_ROOT / "data" / "results" / "figures"
VIZ_DIR.mkdir(parents=True, exist_ok=True)

N_GENOMES = 46


# ── helpers ──────────────────────────────────────────────────────────────────

def load_genome_manifest() -> dict:
    """Return {genome_name: PATHOGEN|COMMENSAL} from the curated manifest."""
    manifest_path = REPO_ROOT / "data" / "genomes" / "genome_manifest.tsv"
    manifest = pd.read_csv(manifest_path, sep="\t")
    manifest = manifest[manifest["label"].notna() & (manifest["label"] != "label")]

    # Binary classification: known pathogenic pathotypes vs non-pathogenic
    pathogenic_types = {"EHEC", "EPEC", "ETEC", "UPEC", "NMEC", "EAEC", "AIEC"}
    result = {}
    for _, row in manifest.iterrows():
        label = str(row["label"]).strip()
        ptype = str(row["pathotype"]).strip().upper()
        result[label] = "PATHOGEN" if ptype in pathogenic_types else "COMMENSAL"
    return result


# ── 1. Pangenome composition ─────────────────────────────────────────────────

def plot_pangenome_composition(gc: pd.DataFrame) -> None:
    # Count how many unique genomes each gene cluster appears in
    cluster_counts = (
        gc.groupby("gene_cluster_id")["genome_name"]
        .nunique()
        .reset_index(name="n_genomes")
    )

    total = len(cluster_counts)

    # Partition into classic pangenome categories
    def category(n):
        if n == N_GENOMES:
            return "Core\n(all 46)"
        elif n >= int(0.95 * N_GENOMES):
            return "Soft-core\n(≥95%)"
        elif n >= 2:
            return "Accessory\n(2–43)"
        else:
            return "Unique\n(1 genome)"

    cluster_counts["category"] = cluster_counts["n_genomes"].apply(category)
    cat_order = ["Core\n(all 46)", "Soft-core\n(≥95%)",
                 "Accessory\n(2–43)", "Unique\n(1 genome)"]
    counts = cluster_counts["category"].value_counts().reindex(cat_order, fill_value=0)

    colours = ["#1565C0", "#42A5F5", "#FFA726", "#EF5350"]

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    # Bar chart
    ax = axes[0]
    bars = ax.bar(counts.index, counts.values, color=colours, edgecolor="white",
                  linewidth=0.8, width=0.6)
    for bar, val in zip(bars, counts.values):
        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 30,
                f"{val:,}\n({val/total*100:.1f}%)",
                ha="center", va="bottom", fontsize=9)
    ax.set_ylabel("Number of gene clusters", fontsize=11)
    ax.set_title(f"Pangenome composition\n46 genomes · {total:,} gene clusters total",
                 fontsize=11, pad=8)
    ax.set_ylim(0, counts.max() * 1.18)
    ax.grid(axis="y", color="#eeeeee", linewidth=0.6)
    ax.spines[["top", "right"]].set_visible(False)

    # Pie chart
    ax2 = axes[1]
    wedge_props = dict(linewidth=1.2, edgecolor="white")
    ax2.pie(counts.values, labels=counts.index, colors=colours,
            autopct="%1.1f%%", startangle=140,
            wedgeprops=wedge_props, textprops={"fontsize": 10})
    ax2.set_title("Gene cluster partition\n(pangenome shell)", fontsize=11, pad=8)

    plt.suptitle("*E. coli* O157:H7 Sakai — 46-genome Anvi'o pangenome",
                 fontsize=12, y=1.01, style="italic")
    plt.tight_layout()
    out = VIZ_DIR / "pangenome_composition.png"
    plt.savefig(out, dpi=180, bbox_inches="tight")
    plt.close()
    log.info("Saved: %s", out)


# ── 2. Functional enrichment dot plot ────────────────────────────────────────

def plot_enrichment(enrich: pd.DataFrame) -> None:
    sig = enrich[enrich["adjusted_q_value"] < 0.05].copy()
    sig = sig[sig["associated_groups"].str.contains("PATHOGEN", na=False)]
    sig["presence_diff"] = sig["p_PATHOGEN"] - sig["p_COMMENSAL"]
    sig = sig.sort_values("enrichment_score", ascending=True)

    # Truncate long names
    def shorten(name, maxlen=52):
        return name if len(name) <= maxlen else name[:maxlen - 1] + "…"

    sig["label"] = sig["COG14_FUNCTION"].apply(shorten)

    fig, ax = plt.subplots(figsize=(9, max(4, len(sig) * 0.65 + 1.2)))

    scatter = ax.scatter(
        sig["presence_diff"],
        range(len(sig)),
        s=sig["enrichment_score"] * 5,
        c=sig["enrichment_score"],
        cmap="YlOrRd",
        vmin=0, vmax=sig["enrichment_score"].max() * 1.1,
        edgecolors="#333333", linewidths=0.5, zorder=3,
    )

    for i, row in enumerate(sig.itertuples()):
        ax.annotate(
            f"q={row.adjusted_q_value:.3f}",
            xy=(row.presence_diff, i),
            xytext=(5, 0), textcoords="offset points",
            va="center", fontsize=7.5, color="#555555",
        )

    ax.set_yticks(range(len(sig)))
    ax.set_yticklabels(sig["label"].tolist(), fontsize=9)
    ax.set_xlabel("Presence difference (pathogenic fraction − commensal fraction)",
                  fontsize=10)
    ax.set_title("PATHOGEN-enriched COG14 functions (q < 0.05)\n"
                 "Anvi'o functional enrichment · 46-genome pangenome",
                 fontsize=11, pad=8)
    ax.axvline(0, color="#aaaaaa", linewidth=0.8, linestyle="--")
    ax.grid(axis="x", color="#eeeeee", linewidth=0.6)
    ax.spines[["top", "right"]].set_visible(False)

    cbar = plt.colorbar(scatter, ax=ax, pad=0.01)
    cbar.set_label("Enrichment score", fontsize=9)

    plt.tight_layout()
    out = VIZ_DIR / "pangenome_enrichment.png"
    plt.savefig(out, dpi=180, bbox_inches="tight")
    plt.close()
    log.info("Saved: %s", out)


# ── 3. Presence/absence heatmap (top 30 differentially present clusters) ─────

def plot_presence_heatmap(gc: pd.DataFrame, enrich: pd.DataFrame,
                          pathotype: dict) -> None:
    # Get significant PATHOGEN-enriched cluster IDs
    sig = enrich[
        (enrich["adjusted_q_value"] < 0.05) &
        (enrich["associated_groups"].str.contains("PATHOGEN", na=False))
    ]
    # Explode gene_clusters_ids column
    sig_clusters = []
    for _, row in sig.iterrows():
        ids = [c.strip() for c in str(row["gene_clusters_ids"]).split(",")]
        for gc_id in ids:
            sig_clusters.append({
                "gene_cluster_id": gc_id,
                "function": row["COG14_FUNCTION"][:40],
            })
    sig_df = pd.DataFrame(sig_clusters).drop_duplicates("gene_cluster_id")

    # Build presence/absence matrix for these clusters
    top_clusters = sig_df["gene_cluster_id"].tolist()[:30]
    sub = gc[gc["gene_cluster_id"].isin(top_clusters)]
    matrix = (
        sub.groupby(["genome_name", "gene_cluster_id"])
        .size()
        .gt(0)
        .unstack(fill_value=False)
        .astype(int)
    )

    # Order genomes: pathogens first
    genomes = list(matrix.index)
    path_genomes = sorted([g for g in genomes if pathotype.get(g) == "PATHOGEN"])
    comm_genomes = sorted([g for g in genomes if pathotype.get(g) != "PATHOGEN"])
    ordered_genomes = path_genomes + comm_genomes
    matrix = matrix.reindex([g for g in ordered_genomes if g in matrix.index])

    if matrix.empty or matrix.shape[1] == 0:
        log.warning("No significant cluster data for heatmap — skipping")
        return

    # Short genome labels
    def short_label(name):
        return name[:20] + "…" if len(name) > 20 else name

    matrix.index = [short_label(g) for g in matrix.index]

    fig, ax = plt.subplots(figsize=(max(8, matrix.shape[1] * 0.55),
                                    max(6, matrix.shape[0] * 0.35)))

    sns.heatmap(
        matrix, ax=ax,
        cmap=["#f0f0f0", "#1565C0"],
        linewidths=0.3, linecolor="#cccccc",
        cbar_kws={"label": "Gene cluster present (1) / absent (0)",
                  "ticks": [0, 1]},
        xticklabels=True, yticklabels=True,
    )

    n_path = len([g for g in matrix.index
                  if any(k.lower() in g.lower()
                         for k in ["O157", "EDL", "EC4", "Sakai", "EPEC",
                                   "CFT", "UTI", "O111", "O103", "O145"])])

    ax.axhline(n_path, color="#E53935", linewidth=1.8, linestyle="--")
    ax.set_title("Presence/absence — PATHOGEN-enriched gene clusters\n"
                 f"Top {matrix.shape[1]} significant clusters · red line = pathogen/commensal boundary",
                 fontsize=11, pad=8)
    ax.set_xlabel("Gene cluster ID", fontsize=10)
    ax.set_ylabel("Genome", fontsize=10)
    ax.tick_params(axis="x", labelsize=7, rotation=45)
    ax.tick_params(axis="y", labelsize=7)

    plt.tight_layout()
    out = VIZ_DIR / "pangenome_heatmap.png"
    plt.savefig(out, dpi=180, bbox_inches="tight")
    plt.close()
    log.info("Saved: %s", out)


def main():
    log.info("Loading gene clusters (%s)…", PAN_DIR / "gene_clusters.tsv")
    gc = pd.read_csv(PAN_DIR / "gene_clusters.tsv", sep="\t")
    log.info("  %d rows, %d unique clusters, %d genomes",
             len(gc),
             gc["gene_cluster_id"].nunique(),
             gc["genome_name"].nunique())

    log.info("Loading enrichment scores…")
    enrich = pd.read_csv(PAN_DIR / "enrichment_scores.tsv", sep="\t")

    pathotype = load_genome_manifest()
    log.info("  %d pathogenic / %d commensal genomes detected",
             sum(v == "PATHOGEN" for v in pathotype.values()),
             sum(v == "COMMENSAL" for v in pathotype.values()))

    log.info("Plotting pangenome composition…")
    plot_pangenome_composition(gc)

    log.info("Plotting functional enrichment…")
    plot_enrichment(enrich)

    log.info("Plotting presence/absence heatmap…")
    plot_presence_heatmap(gc, enrich, pathotype)

    log.info("All pangenome figures saved to %s", VIZ_DIR)


if __name__ == "__main__":
    main()
