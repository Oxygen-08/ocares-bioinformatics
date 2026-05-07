#!/usr/bin/env python3
"""
Plot alignment block distribution across the O157:H7 Sakai reference
for all 45 comparison strains, coloured by tier.

Each horizontal segment represents one NUCmer alignment block.
Strains are ordered by pathotype group, then alphabetically.
"""

import csv
import logging
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import pandas as pd

logging.basicConfig(level=logging.INFO, format="%(asctime)s  %(message)s", datefmt="%H:%M:%S")
log = logging.getLogger(__name__)

REPO_ROOT   = Path(__file__).parents[2]
RESULTS_DIR = REPO_ROOT / "data" / "results"
VIZ_DIR     = RESULTS_DIR / "figures"
VIZ_DIR.mkdir(parents=True, exist_ok=True)

TIER_COLORS = {
    "CONSERVED": "#2196F3",   # blue
    "MODERATE":  "#FF9800",   # orange
    "DIVERGED":  "#F44336",   # red
}

# Pathotype display order and group colours for y-axis labels
PATHOTYPE_ORDER = ["EHEC", "UPEC", "ETEC", "EAEC", "EPEC", "NMEC", "AIEC", "K12", "commensal"]
PATHOTYPE_LABEL_COLOR = {
    "EHEC": "#b71c1c", "UPEC": "#880e4f", "ETEC": "#4a148c",
    "EAEC": "#1a237e", "EPEC": "#01579b", "NMEC": "#006064",
    "AIEC": "#1b5e20", "K12":  "#33691e", "commensal": "#827717",
}


def load_blocks() -> pd.DataFrame:
    path = RESULTS_DIR / "tiered_blocks.tsv"
    df = pd.read_csv(path, sep="\t")
    return df


def strain_sort_key(row):
    order = PATHOTYPE_ORDER.index(row["pathotype"]) if row["pathotype"] in PATHOTYPE_ORDER else 99
    return (order, row["label"])


def main() -> None:
    df = load_blocks()
    log.info("Loaded %d alignment blocks from %d strains", len(df), df["label"].nunique())

    # Build ordered strain list
    strain_meta = df[["label", "pathotype"]].drop_duplicates()
    strain_meta = strain_meta.sort_values(
        by=["pathotype", "label"],
        key=lambda col: col.map(
            lambda v: PATHOTYPE_ORDER.index(v) if v in PATHOTYPE_ORDER else 99
        ) if col.name == "pathotype" else col
    )
    strains = strain_meta["label"].tolist()
    y_pos   = {s: i for i, s in enumerate(strains)}
    n       = len(strains)

    # Reference genome length from max ref_end
    ref_len = df["ref_end"].max()
    log.info("Reference length (max ref_end): %d bp", ref_len)

    fig, ax = plt.subplots(figsize=(18, max(10, n * 0.28)))

    for _, row in df.iterrows():
        y     = y_pos[row["label"]]
        color = TIER_COLORS.get(row["tier"], "#999999")
        alpha = 0.55 if row["tier"] == "CONSERVED" else 0.80
        lw    = 0.7 if row["tier"] == "CONSERVED" else 1.1
        ax.plot(
            [row["ref_start"], row["ref_end"]],
            [y, y],
            color=color, linewidth=lw, alpha=alpha, solid_capstyle="butt"
        )

    # Y-axis labels coloured by pathotype
    ax.set_yticks(range(n))
    labels = []
    for s in strains:
        pt = strain_meta.loc[strain_meta["label"] == s, "pathotype"].values[0]
        labels.append(s)
    ax.set_yticklabels(labels, fontsize=7)

    # Colour y-tick labels by pathotype
    for tick, s in zip(ax.get_yticklabels(), strains):
        pt = strain_meta.loc[strain_meta["label"] == s, "pathotype"].values[0]
        tick.set_color(PATHOTYPE_LABEL_COLOR.get(pt, "#333333"))

    # X-axis: Mb
    ax.set_xlabel("O157:H7 Sakai reference position (Mb)", fontsize=11)
    ax.set_xlim(0, ref_len)
    ax.xaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: f"{x/1e6:.1f}"))
    ax.set_ylabel("Strain", fontsize=11)
    ax.set_title(
        "Alignment block distribution across O157:H7 Sakai reference\n"
        "45 comparison strains · coloured by divergence tier",
        fontsize=13, pad=12,
    )

    # Tier legend
    handles = [
        mpatches.Patch(color=TIER_COLORS["CONSERVED"], label="Tier I — Conserved (≥95% identity)"),
        mpatches.Patch(color=TIER_COLORS["MODERATE"],  label="Tier II — Moderate (85–95% identity)"),
        mpatches.Patch(color=TIER_COLORS["DIVERGED"],  label="Tier III — Diverged (<85% identity)"),
    ]
    ax.legend(handles=handles, loc="upper right", fontsize=9, framealpha=0.9)

    # Pathotype group separators
    current_pt = None
    for i, s in enumerate(strains):
        pt = strain_meta.loc[strain_meta["label"] == s, "pathotype"].values[0]
        if pt != current_pt:
            if i > 0:
                ax.axhline(i - 0.5, color="#cccccc", linewidth=0.6, linestyle="--")
            current_pt = pt
            ax.text(ref_len * 1.005, i, pt, fontsize=7, va="center",
                    color=PATHOTYPE_LABEL_COLOR.get(pt, "#333333"), fontweight="bold")

    ax.invert_yaxis()
    ax.grid(axis="x", color="#eeeeee", linewidth=0.5)
    plt.tight_layout()

    out = VIZ_DIR / "alignment_landscape.png"
    plt.savefig(out, dpi=180, bbox_inches="tight")
    plt.close()
    log.info("Saved: %s", out)


if __name__ == "__main__":
    main()
