#!/usr/bin/env python3
"""
Biological sanity check: validate that DIVERGED markers map to known
O157:H7 Sakai pathogenicity islands and accessory genome regions.

Two independent tests:
  1. K12 absence test — for each DIVERGED marker, check if K12 strains
     (commensal lab strains) have NUCmer alignment coverage.  Absence = truly
     accessory genome (PAI / O-island).
  2. COG functional enrichment overlap — do DIVERGED markers sit in gene
     clusters carrying pathogen-enriched COG functions?

Outputs:
  data/results/bio_sanity/diverged_marker_summary.tsv
  data/results/bio_sanity/k12_absence_rate.txt
  data/results/figures/bio_sanity_overview.png
"""

import logging
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd

logging.basicConfig(level=logging.INFO, format="%(asctime)s  %(message)s",
                    datefmt="%H:%M:%S")
log = logging.getLogger(__name__)

REPO_ROOT   = Path(__file__).parents[2]
RESULTS_DIR = REPO_ROOT / "data" / "results"
OUT_DIR     = RESULTS_DIR / "bio_sanity"
VIZ_DIR     = RESULTS_DIR / "figures"
OUT_DIR.mkdir(parents=True, exist_ok=True)
VIZ_DIR.mkdir(parents=True, exist_ok=True)

# K12 strains = commensal lab strains (no known PAIs)
K12_LABELS = {"K12_MG1655", "K12_DH10B", "K12_W3110", "BL21", "BL21_DE3",
              "B_REL606", "MG1655", "W"}

# Approximate O-island / virulence locus coordinates on NC_002695.2
# Hayashi et al. 2001 (PMID 11085612), Ohnishi et al. 2002
KNOWN_LOCI = [
    {"name": "Stx1 prophage (SpLE2)",   "start":  1_378_635, "end":  1_450_000},
    {"name": "Stx2 prophage (SpLE3)",   "start":  1_450_000, "end":  1_540_000},
    {"name": "Tellurite resistance",     "start":  1_189_000, "end":  1_235_000},
    {"name": "OI-43 (type III related)", "start":  2_458_000, "end":  2_545_000},
    {"name": "OI-48 (TTSS-2)",          "start":  2_750_000, "end":  2_820_000},
    {"name": "LEE (OI-148)",            "start":  5_030_000, "end":  5_100_000},
    {"name": "Flagellar region",         "start":  1_940_000, "end":  2_010_000},
]


def load_data():
    meta   = pd.read_csv(RESULTS_DIR / "markers" / "marker_metadata.tsv", sep="\t")
    blocks = pd.read_csv(RESULTS_DIR / "tiered_blocks.tsv", sep="\t")
    enrich = pd.read_csv(RESULTS_DIR / "pangenome" / "enrichment_scores.tsv", sep="\t")
    return meta, blocks, enrich


def k12_absence_test(meta: pd.DataFrame, blocks: pd.DataFrame) -> pd.DataFrame:
    """For each DIVERGED marker, is it absent from all K12 strains?"""
    diverged = meta[meta["tier"] == "DIVERGED"].copy()
    k12_blocks = blocks[blocks["label"].isin(K12_LABELS)]

    rows = []
    for _, m in diverged.iterrows():
        overlap = k12_blocks[
            (k12_blocks["ref_contig"] == m["contig"]) &
            (k12_blocks["ref_start"]  <= m["end"])    &
            (k12_blocks["ref_end"]    >= m["start"])
        ]
        k12_strains_present = overlap["label"].nunique()
        rows.append({
            "marker_id":          m["marker_id"],
            "contig":             m["contig"],
            "start":              m["start"],
            "end":                m["end"],
            "length":             m["length"],
            "k12_strains_with_coverage": k12_strains_present,
            "absent_from_k12":    k12_strains_present == 0,
        })
    return pd.DataFrame(rows)


def known_locus_overlap(df: pd.DataFrame) -> pd.DataFrame:
    """Flag markers that overlap known virulence loci."""
    results = []
    for locus in KNOWN_LOCI:
        hits = df[
            (df["contig"] == "NC_002695.2") &
            (df["start"]  <= locus["end"])  &
            (df["end"]    >= locus["start"])
        ]["marker_id"].tolist()
        if hits:
            results.append({"locus": locus["name"], "n_markers": len(hits),
                             "marker_ids": ", ".join(hits[:5])})
    return pd.DataFrame(results)


def enrichment_overlap(df: pd.DataFrame, enrich: pd.DataFrame) -> int:
    """Count DIVERGED markers whose gene clusters carry enriched COG functions."""
    gc_path = RESULTS_DIR / "pangenome" / "gene_clusters.tsv"
    gene_calls_path = RESULTS_DIR / "pangenome" / "sakai_gene_calls.tsv"
    if not gc_path.exists() or not gene_calls_path.exists():
        return 0

    # Build set of enriched gene cluster IDs (q < 0.05, PATHOGEN)
    sig = enrich[(enrich["adjusted_q_value"] < 0.05) &
                 (enrich["associated_groups"] == "PATHOGEN")]
    enriched_clusters = set()
    for _, row in sig.iterrows():
        for cid in str(row["gene_clusters_ids"]).split(","):
            enriched_clusters.add(cid.strip())

    from Bio import SeqIO
    ref_fasta     = REPO_ROOT / "data" / "genomes" / "O157H7_Sakai" / "O157H7_Sakai.fna"
    reform_fasta  = RESULTS_DIR / "pangenome" / "reformatted" / "O157H7_Sakai.fna"
    orig_ids      = [r.id for r in SeqIO.parse(str(ref_fasta), "fasta")]
    reformed_ids  = [r.id for r in SeqIO.parse(str(reform_fasta), "fasta")]
    r2o = dict(zip(reformed_ids, orig_ids))

    gc    = pd.read_csv(gc_path, sep="\t",
                        usecols=["gene_caller_id", "gene_cluster_id", "genome_name"])
    sakai = gc[gc["genome_name"] == "O157H7_Sakai"][["gene_caller_id", "gene_cluster_id"]]
    calls = pd.read_csv(gene_calls_path, sep="\t",
                        usecols=["gene_callers_id", "contig", "start", "stop"])
    calls["orig_contig"] = calls["contig"].map(r2o)
    calls = calls.merge(sakai, left_on="gene_callers_id",
                        right_on="gene_caller_id", how="left")

    count = 0
    for _, m in df.iterrows():
        overlap = calls[
            (calls["orig_contig"] == m["contig"]) &
            (calls["start"]       <= m["end"])    &
            (calls["stop"]        >= m["start"])  &
            calls["gene_cluster_id"].notna()
        ]
        if any(cid in enriched_clusters for cid in overlap["gene_cluster_id"]):
            count += 1
    return count


def plot_overview(df: pd.DataFrame, locus_df: pd.DataFrame) -> None:
    ref_len  = 5_498_578  # NC_002695.2 length
    fig, axes = plt.subplots(2, 1, figsize=(16, 8),
                             gridspec_kw={"height_ratios": [3, 1]})

    ax = axes[0]
    absent  = df[df["absent_from_k12"]]
    present = df[~df["absent_from_k12"]]
    ax.scatter(absent["start"],  np.ones(len(absent))  * 1,
               c="#F44336", s=20, alpha=0.7, label=f"Absent from K12 (n={len(absent)})")
    ax.scatter(present["start"], np.ones(len(present)) * 0,
               c="#2196F3", s=20, alpha=0.5, label=f"Present in ≥1 K12 (n={len(present)})")

    for locus in KNOWN_LOCI:
        ax.axvspan(locus["start"], locus["end"], alpha=0.12, color="#FF9800")
        ax.text((locus["start"] + locus["end"]) / 2, 1.35,
                locus["name"].split("(")[0].strip(),
                ha="center", va="bottom", fontsize=6.5, color="#E65100", rotation=30)

    ax.set_xlim(0, ref_len)
    ax.set_yticks([0, 1])
    ax.set_yticklabels(["In K12", "Absent from K12"], fontsize=9)
    ax.xaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: f"{x/1e6:.1f}"))
    ax.set_xlabel("O157:H7 Sakai position (Mb)", fontsize=10)
    ax.set_title("DIVERGED marker distribution — K12 absence test\n"
                 "Orange shading = known virulence loci (Hayashi et al. 2001)",
                 fontsize=12, pad=10)
    ax.legend(fontsize=9, loc="upper left")
    ax.grid(axis="x", color="#eeeeee", linewidth=0.5)

    ax2 = axes[1]
    if not locus_df.empty:
        ax2.barh(locus_df["locus"], locus_df["n_markers"], color="#7B1FA2", alpha=0.8)
        ax2.set_xlabel("DIVERGED markers overlapping locus", fontsize=9)
        ax2.set_title("Known virulence loci with DIVERGED marker overlap", fontsize=10)
        ax2.grid(axis="x", color="#eeeeee", linewidth=0.5)
    else:
        ax2.text(0.5, 0.5, "No overlap with known loci",
                 ha="center", va="center", transform=ax2.transAxes, fontsize=10)

    plt.tight_layout()
    out = VIZ_DIR / "bio_sanity_overview.png"
    plt.savefig(out, dpi=180, bbox_inches="tight")
    plt.close()
    log.info("Saved: %s", out)


def main():
    meta, blocks, enrich = load_data()
    log.info("Loaded %d markers (%d DIVERGED)",
             len(meta), (meta["tier"] == "DIVERGED").sum())

    log.info("Running K12 absence test...")
    absence_df = k12_absence_test(meta, blocks)
    n_absent   = absence_df["absent_from_k12"].sum()
    pct_absent = 100 * n_absent / len(absence_df)
    log.info("  %d / %d DIVERGED markers (%.1f%%) absent from all K12 strains",
             n_absent, len(absence_df), pct_absent)

    log.info("Checking known virulence locus overlap...")
    locus_df = known_locus_overlap(absence_df)
    if not locus_df.empty:
        log.info("\n%s", locus_df.to_string(index=False))

    log.info("Checking COG-enriched cluster overlap...")
    n_enrich = enrichment_overlap(absence_df, enrich)
    log.info("  %d DIVERGED markers overlap COG-enriched gene clusters", n_enrich)

    # Save summary
    absence_df.to_csv(OUT_DIR / "diverged_marker_summary.tsv", sep="\t", index=False)
    summary_txt = (
        f"DIVERGED markers: {len(absence_df)}\n"
        f"Absent from all K12 strains: {n_absent} ({pct_absent:.1f}%)\n"
        f"Overlap known virulence loci: {locus_df['n_markers'].sum() if not locus_df.empty else 0}\n"
        f"Overlap COG-enriched clusters: {n_enrich}\n"
    )
    (OUT_DIR / "k12_absence_rate.txt").write_text(summary_txt)
    log.info("\n%s", summary_txt)

    plot_overview(absence_df, locus_df)
    log.info("Bio sanity check complete. Results in %s", OUT_DIR)


if __name__ == "__main__":
    main()
