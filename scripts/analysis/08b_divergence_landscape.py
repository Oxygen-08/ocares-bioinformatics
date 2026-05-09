"""
Genome-coordinate divergence landscape plots for the minimap2-gradient branch.

Produces three figures:
  1. divergence_landscape_chromosome.png
       Full Sakai chromosome (NC_002695.2) coloured by mean divergence score.
       HIGH-gradient candidate regions shaded. Known PAI landmarks annotated.

  2. divergence_landscape_plasmids.png
       pO157 megaplasmid (NC_002128.1) and F-plasmid (NC_002127.1) divergence.

  3. divergence_confidence_tiers.png
       Chromosome with NUCmer markers coloured by confidence tier:
         Tier 1 (green)  — both NUCmer + minimap2 agree (333 markers)
         Tier 2 (orange) — NUCmer-only (82 markers)
       Positive class (DIVERGED) shown as filled, MODERATE as outline.
"""

from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

REPO_ROOT  = Path(__file__).resolve().parents[2]
MM_DIR     = REPO_ROOT / "data" / "results" / "minimap2"
ML_DIR     = REPO_ROOT / "data" / "results" / "ml"
FIG_DIR    = REPO_ROOT / "data" / "results" / "figures"
FIG_DIR.mkdir(parents=True, exist_ok=True)

# Known Sakai PAI coordinates (NC_002695.2) from Hayashi et al. 2001 PMID 11258796
# and Perna et al. 2001 PMID 11206551 — approximate positions for annotation only
KNOWN_PAIS = [
    {"name": "LEE",       "start": 2_796_000, "end": 2_852_000, "color": "#c0392b"},
    {"name": "OI-1",      "start":    17_000, "end":    68_000, "color": "#8e44ad"},
    {"name": "OI-7",      "start":   516_000, "end":   562_000, "color": "#8e44ad"},
    {"name": "OI-36",     "start": 1_640_000, "end": 1_693_000, "color": "#8e44ad"},
    {"name": "OI-43/48",  "start": 2_100_000, "end": 2_195_000, "color": "#d35400"},
    {"name": "OI-57",     "start": 2_450_000, "end": 2_510_000, "color": "#8e44ad"},
    {"name": "Sp1",       "start": 1_310_000, "end": 1_360_000, "color": "#27ae60"},
    {"name": "Sp4",       "start": 2_560_000, "end": 2_610_000, "color": "#27ae60"},
    {"name": "Sp11",      "start": 4_500_000, "end": 4_560_000, "color": "#27ae60"},
    {"name": "Sp15",      "start": 5_100_000, "end": 5_160_000, "color": "#27ae60"},
]

SCHEME_COLORS = {
    "LOW":  "#3498db",
    "MID":  "#f39c12",
    "HIGH": "#e74c3c",
}


def classify_score(score, low_hi=0.20, mid_hi=0.60):
    if score <= low_hi:
        return "LOW"
    elif score <= mid_hi:
        return "MID"
    return "HIGH"


def load_data():
    windows = pd.read_csv(MM_DIR / "window_scores_500bp.tsv", sep="\t")
    windows["gradient"] = windows["mean_div_score"].apply(classify_score)
    candidates = pd.read_csv(MM_DIR / "candidates_schemeA_500bp.tsv", sep="\t")
    markers = pd.read_csv(ML_DIR / "feature_matrix.tsv", sep="\t")
    return windows, candidates, markers


# ── Figure 1: Chromosome divergence landscape ──────────────────────────────────

def plot_chromosome_landscape(windows, candidates):
    chrom = windows[windows["contig"] == "NC_002695.2"].copy()
    chrom_cands = candidates[candidates["contig"] == "NC_002695.2"].copy()
    chrom_len = chrom["w_end"].max()

    fig, axes = plt.subplots(3, 1, figsize=(16, 9),
                             gridspec_kw={"height_ratios": [3, 0.6, 0.5]})
    fig.suptitle("Sakai Chromosome (NC_002695.2) — minimap2 Divergence Landscape\n"
                 "30-commensal panel, Scheme A (LOW ≤0.20, MID ≤0.60, HIGH >0.60)",
                 fontsize=11, y=0.98)

    # ── Axes 0: divergence score line + gradient fill ──
    ax = axes[0]
    x_mid = (chrom["w_start"] + chrom["w_end"]) / 2 / 1e6

    ax.fill_between(x_mid, 0, chrom["mean_div_score"],
                    where=chrom["gradient"] == "LOW",
                    color=SCHEME_COLORS["LOW"], alpha=0.6, label="LOW (≤0.20)")
    ax.fill_between(x_mid, 0, chrom["mean_div_score"],
                    where=chrom["gradient"] == "MID",
                    color=SCHEME_COLORS["MID"], alpha=0.6, label="MID (0.20–0.60)")
    ax.fill_between(x_mid, 0, chrom["mean_div_score"],
                    where=chrom["gradient"] == "HIGH",
                    color=SCHEME_COLORS["HIGH"], alpha=0.8, label="HIGH (>0.60)")

    ax.axhline(0.20, color="#3498db", lw=0.8, ls="--", alpha=0.5)
    ax.axhline(0.60, color="#e74c3c", lw=0.8, ls="--", alpha=0.5)
    ax.set_ylim(0, 1.05)
    ax.set_ylabel("Mean divergence score\n(1 − coverage × identity)", fontsize=9)
    ax.set_xlim(0, chrom_len / 1e6)
    ax.legend(loc="upper right", fontsize=8, framealpha=0.8)
    ax.tick_params(labelbottom=False)

    # Annotate known PAIs
    for pai in KNOWN_PAIS:
        mid = (pai["start"] + pai["end"]) / 2 / 1e6
        ax.axvspan(pai["start"] / 1e6, pai["end"] / 1e6,
                   color=pai["color"], alpha=0.15, zorder=0)
        ax.text(mid, 1.02, pai["name"], ha="center", va="bottom",
                fontsize=6, color=pai["color"], rotation=45)

    # ── Axes 1: candidate regions ──
    ax1 = axes[1]
    ax1.set_xlim(0, chrom_len / 1e6)
    ax1.set_ylim(0, 1)
    ax1.set_yticks([])
    ax1.set_ylabel("Candidates\n(Scheme A)", fontsize=8)
    ax1.tick_params(labelbottom=False)

    for _, row in chrom_cands.iterrows():
        ax1.axvspan(row["region_start"] / 1e6, row["region_end"] / 1e6,
                    color="#e74c3c", alpha=0.7)

    # ── Axes 2: negative control comparison ──
    ax2 = axes[2]
    ax2.set_xlim(0, chrom_len / 1e6)
    neg_dir = MM_DIR / "negative_control"
    neg_paf_scores = []
    if neg_dir.exists():
        neg_tsv = neg_dir / "neg_control_window_scores.tsv"
        if neg_tsv.exists():
            neg_w = pd.read_csv(neg_tsv, sep="\t")
            neg_chrom = neg_w[neg_w["contig"] == "NC_002695.2"]
            x_neg = (neg_chrom["w_start"] + neg_chrom["w_end"]) / 2 / 1e6
            ax2.fill_between(x_neg, 0, neg_chrom["mean_div_score"],
                             color="#7f8c8d", alpha=0.6, label="Pathogen vs pathogen")
            ax2.set_ylim(0, 1.05)
            ax2.legend(loc="upper right", fontsize=7)

    ax2.set_ylabel("Neg. control\ndivergence", fontsize=8)
    ax2.set_xlabel("Chromosomal position (Mb)", fontsize=9)
    ax2.set_ylim(0, 1.05)

    plt.tight_layout(rect=[0, 0, 1, 0.97])
    out = FIG_DIR / "divergence_landscape_chromosome.png"
    plt.savefig(out, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"Saved: {out}")


# ── Figure 2: Plasmid landscapes ──────────────────────────────────────────────

def plot_plasmid_landscapes(windows):
    plasmids = {
        "NC_002128.1": ("pO157 megaplasmid", 92721),
        "NC_002127.1": ("F-plasmid (pOSAK1)", 3306),
    }

    fig, axes = plt.subplots(len(plasmids), 1, figsize=(14, 5))
    fig.suptitle("Sakai Plasmids — minimap2 Divergence Landscape", fontsize=11)

    for ax, (contig, (name, length)) in zip(axes, plasmids.items()):
        plas = windows[windows["contig"] == contig].copy()
        if plas.empty:
            ax.text(0.5, 0.5, f"{name}: no data", ha="center", transform=ax.transAxes)
            continue
        x_mid = (plas["w_start"] + plas["w_end"]) / 2 / 1e3

        for state, color in SCHEME_COLORS.items():
            mask = plas["gradient"] == state
            ax.fill_between(x_mid, 0, plas["mean_div_score"].where(mask, 0),
                            color=color, alpha=0.7, label=state)

        ax.axhline(0.20, color="#3498db", lw=0.8, ls="--", alpha=0.5)
        ax.axhline(0.60, color="#e74c3c", lw=0.8, ls="--", alpha=0.5)
        ax.set_ylim(0, 1.05)
        ax.set_ylabel("Mean divergence", fontsize=8)
        ax.set_xlabel("Position (kb)", fontsize=8)
        ax.set_title(f"{name} ({contig})", fontsize=9)
        ax.legend(loc="upper right", fontsize=7)

    plt.tight_layout()
    out = FIG_DIR / "divergence_landscape_plasmids.png"
    plt.savefig(out, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"Saved: {out}")


# ── Figure 3: Confidence tier map ─────────────────────────────────────────────

def plot_confidence_tiers(windows, markers):
    chrom = windows[windows["contig"] == "NC_002695.2"].copy()
    chrom_len = chrom["w_end"].max()
    chrom_markers = markers[markers["contig"] == "NC_002695.2"].copy()

    # Tier assignment
    chrom_markers = chrom_markers.copy()
    chrom_markers["tier1"] = chrom_markers["mean_divergence"] > 0

    fig, axes = plt.subplots(2, 1, figsize=(16, 6),
                             gridspec_kw={"height_ratios": [2, 1]})
    fig.suptitle("Sakai Chromosome — NUCmer Marker Confidence Tiers\n"
                 "Tier 1: both NUCmer + minimap2 agree (333) | Tier 2: NUCmer-only (82)",
                 fontsize=10)

    # Background divergence
    ax = axes[0]
    x_mid = (chrom["w_start"] + chrom["w_end"]) / 2 / 1e6
    ax.fill_between(x_mid, 0, chrom["mean_div_score"],
                    color="#ecf0f1", alpha=0.9)
    ax.fill_between(x_mid, 0, chrom["mean_div_score"],
                    where=chrom["gradient"] == "HIGH",
                    color="#fadbd8", alpha=0.8)

    # Plot markers by tier and class
    for _, row in chrom_markers.iterrows():
        mid   = (row["start"] + row["end"]) / 2 / 1e6
        score = row.get("mean_divergence", 0)
        is_div = row["tier"] == "DIVERGED"
        color  = "#27ae60" if row["tier1"] else "#e67e22"
        marker = "D" if is_div else "o"
        size   = 40 if is_div else 20
        ax.scatter(mid, score + 0.03, c=color, marker=marker,
                   s=size, zorder=5, edgecolors="white", linewidths=0.5)

    ax.set_xlim(0, chrom_len / 1e6)
    ax.set_ylim(0, 1.15)
    ax.set_ylabel("Mean divergence score", fontsize=9)
    ax.tick_params(labelbottom=False)

    legend_elements = [
        mpatches.Patch(color="#27ae60", label="Tier 1 — both pipelines agree"),
        mpatches.Patch(color="#e67e22", label="Tier 2 — NUCmer only"),
        plt.scatter([], [], marker="D", c="grey", s=40, label="DIVERGED (positive class)"),
        plt.scatter([], [], marker="o", c="grey", s=20, label="MODERATE (negative class)"),
    ]
    ax.legend(handles=legend_elements, loc="upper right", fontsize=8, framealpha=0.9)

    # Tier bar
    ax2 = axes[1]
    ax2.set_xlim(0, chrom_len / 1e6)
    ax2.set_ylim(0, 2)
    ax2.set_yticks([])
    ax2.set_ylabel("Tier", fontsize=8)
    ax2.set_xlabel("Chromosomal position (Mb)", fontsize=9)

    for _, row in chrom_markers.iterrows():
        s = row["start"] / 1e6
        e = row["end"] / 1e6
        color = "#27ae60" if row["tier1"] else "#e67e22"
        y = 1.0 if row["tier"] == "DIVERGED" else 0.2
        h = 0.7 if row["tier"] == "DIVERGED" else 0.5
        ax2.barh(y, e - s, left=s, height=h, color=color, alpha=0.8)

    plt.tight_layout(rect=[0, 0, 1, 0.96])
    out = FIG_DIR / "divergence_confidence_tiers.png"
    plt.savefig(out, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"Saved: {out}")


if __name__ == "__main__":
    print("Loading data...")
    windows, candidates, markers = load_data()

    print("Plotting chromosome landscape...")
    plot_chromosome_landscape(windows, candidates)

    print("Plotting plasmid landscapes...")
    plot_plasmid_landscapes(windows)

    print("Plotting confidence tiers...")
    plot_confidence_tiers(windows, markers)

    print("Done. All figures saved to:", FIG_DIR)
