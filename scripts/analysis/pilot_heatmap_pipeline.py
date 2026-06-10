#!/usr/bin/env python3
"""
pilot_heatmap_pipeline.py
─────────────────────────
Replicates the original pilot study:
  1. NUCmer Sakai vs SE11 → divergent (Sakai-specific) marker regions
  2. BLAST those markers against a 23-strain panel
  3. Call presence / absence per strain
  4. Render a publication-quality heatmap

Panel composition
  - 10 pathogenic E. coli  (from thesis manifest)
  - 10 commensal / lab E. coli  (from thesis manifest)
  - E. fergusonii EFCF056  (other Escherichia species)
  - E. albertii 6S-65-1    (other Escherichia species)
  - S. aureus N315         (other genus / class outgroup)
"""

import os, subprocess, sys, textwrap
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.gridspec as gridspec
from matplotlib.colors import ListedColormap
from matplotlib.ticker import MaxNLocator

# ── directories ──────────────────────────────────────────────────────────────
GENOME_DIR  = "/Users/user/bioinformatics/data/genomes"
RAW_GENOMES = "/Users/user/Documents/MISSION/Data_analytics/Bioinformatics/ucd_project/ucd_project2/data/raw_genomes"
WORK_DIR    = "/Users/user/bioinformatics/data/results/pilot_heatmap"
FIG_DIR     = "/Users/user/bioinformatics/data/results/figures"
os.makedirs(WORK_DIR, exist_ok=True)
os.makedirs(FIG_DIR,  exist_ok=True)

SAKAI = f"{GENOME_DIR}/O157H7_Sakai/O157H7_Sakai.fna"
SE11  = f"{GENOME_DIR}/SE11/SE11.fna"

# ── panel definition ─────────────────────────────────────────────────────────
PATHOGENS = {
    "O157H7_EDL933":  ("EHEC O157:H7 EDL933",  "EHEC",  f"{GENOME_DIR}/O157H7_EDL933/O157H7_EDL933.fna"),
    "O157H7_EC4115":  ("EHEC O157:H7 EC4115",   "EHEC",  f"{GENOME_DIR}/O157H7_EC4115/O157H7_EC4115.fna"),
    "O157H7_TW14359": ("EHEC O157:H7 TW14359",  "EHEC",  f"{GENOME_DIR}/O157H7_TW14359/O157H7_TW14359.fna"),
    "HUSEC2011":      ("EHEC HUSEC2011",         "EHEC",  f"{GENOME_DIR}/HUSEC2011/HUSEC2011.fna"),
    "O26H11_11368":   ("EHEC O26:H11",           "EHEC",  f"{GENOME_DIR}/O26H11_11368/O26H11_11368.fna"),
    "O103H2_12009":   ("EHEC O103:H2",           "EHEC",  f"{GENOME_DIR}/O103H2_12009/O103H2_12009.fna"),
    "CFT073":         ("UPEC CFT073",            "UPEC",  f"{GENOME_DIR}/CFT073/CFT073.fna"),
    "H10407":         ("ETEC H10407",            "ETEC",  f"{GENOME_DIR}/H10407/H10407.fna"),
    "E2348_69":       ("EPEC E2348/69",          "EPEC",  f"{GENOME_DIR}/E2348_69/E2348_69.fna"),
    "APEC_O1":        ("APEC O1",                "APEC",  f"{GENOME_DIR}/APEC_O1/APEC_O1.fna"),
}

COMMENSALS = {
    "SE11":       ("Commensal SE11",        "COMMENSAL", f"{GENOME_DIR}/SE11/SE11.fna"),
    "HS":         ("Commensal HS",          "COMMENSAL", f"{GENOME_DIR}/HS/HS.fna"),
    "IAI1":       ("Commensal IAI1",        "COMMENSAL", f"{GENOME_DIR}/IAI1/IAI1.fna"),
    "ED1a":       ("Commensal ED1a",        "COMMENSAL", f"{GENOME_DIR}/ED1a/ED1a.fna"),
    "SE15":       ("Commensal SE15",        "COMMENSAL", f"{GENOME_DIR}/SE15/SE15.fna"),
    "K12_MG1655": ("K-12 MG1655",          "K12",       f"{GENOME_DIR}/K12_MG1655/K12_MG1655.fna"),
    "K12_W3110":  ("K-12 W3110",           "K12",       f"{GENOME_DIR}/K12_W3110/K12_W3110.fna"),
    "BL21":       ("Lab-B BL21",           "LAB_B",     f"{GENOME_DIR}/BL21/BL21.fna"),
    "Nissle1917": ("Probiotic Nissle 1917","PROBIOTIC",  f"{GENOME_DIR}/Nissle1917/Nissle1917.fna"),
    "ATCC25922":  ("Commensal ATCC25922",  "COMMENSAL", f"{GENOME_DIR}/ATCC25922/ATCC25922.fna"),
}

OUTGROUPS = {
    "E_fergusonii": ("E. fergusonii EFCF056", "OTHER_ESCH", f"{RAW_GENOMES}/efergusonii_efcf056.fasta"),
    "E_albertii":   ("E. albertii 6S-65-1",   "OTHER_ESCH", f"{RAW_GENOMES}/ealberti_6s651.fasta"),
    "S_aureus":     ("S. aureus N315",         "OTHER_GENUS",f"{RAW_GENOMES}/staph_aureus.fasta"),
}

ALL_STRAINS  = {**PATHOGENS, **COMMENSALS, **OUTGROUPS}
STRAIN_ORDER = list(PATHOGENS) + list(COMMENSALS) + list(OUTGROUPS)

def run(cmd, desc=""):
    print(f"  ▶ {desc or cmd[:80]}")
    r = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    if r.returncode != 0:
        print(f"  [STDERR] {r.stderr[:400]}")
        raise RuntimeError(f"Command failed: {cmd[:80]}")
    return r.stdout

# ════════════════════════════════════════════════════════════════════════════
# STEP 1 — NUCmer Sakai vs SE11
# ════════════════════════════════════════════════════════════════════════════
print("\n[1] NUCmer: Sakai vs SE11")
prefix = f"{WORK_DIR}/sakai_vs_se11"
delta  = f"{prefix}.delta"
fdelta = f"{prefix}_filtered.delta"
coords = f"{prefix}_filtered.coords"
aln_bed= f"{WORK_DIR}/sakai_aligned.bed"
div_bed= f"{WORK_DIR}/sakai_divergent.bed"
div_fna= f"{WORK_DIR}/sakai_divergent.fna"
div_fil= f"{WORK_DIR}/sakai_divergent_filtered.fna"
fai    = f"{WORK_DIR}/sakai.fna.fai"
gsizes = f"{WORK_DIR}/sakai_genome.sizes"

if not os.path.exists(delta):
    run(f"nucmer --prefix={prefix} {SAKAI} {SE11}", "nucmer alignment")
else:
    print("  (nucmer delta exists, skipping)")

run(f"delta-filter -m -o 0 {delta} > {fdelta}", "delta-filter (all identity, 1-to-1 mapping)")
run(f"show-coords -rclT {fdelta} > {coords}", "show-coords")

# Genome size index
import shutil
shutil.copy(SAKAI, f"{WORK_DIR}/sakai.fna")
run(f"samtools faidx {WORK_DIR}/sakai.fna", "samtools faidx")
run(f"awk '{{print $1\"\\t\"$2}}' {fai} > {gsizes}", "genome sizes")

# BED of aligned regions in Sakai (col 12 = ref_name, col 1 = ref_start, col 2 = ref_end)
run(
    f"awk 'BEGIN{{OFS=\"\\t\"}} /^[0-9]/{{print $12, $1-1, $2}}' {coords} "
    f"| bedtools sort -i - -g {gsizes} > {aln_bed}",
    "aligned BED"
)

# Complement = divergent (Sakai-specific) regions
run(
    f"bedtools merge -i {aln_bed} | "
    f"bedtools complement -i - -g {gsizes} > {div_bed}",
    "bedtools complement → divergent regions"
)

# Count divergent regions
with open(div_bed) as fh:
    div_regions = [l.strip() for l in fh if l.strip()]
print(f"  ✓ {len(div_regions)} divergent regions in Sakai vs SE11")

# Extract FASTA of divergent regions
run(f"bedtools getfasta -fi {WORK_DIR}/sakai.fna -bed {div_bed} -fo {div_fna}", "getfasta")

# Filter short regions (< 100 bp)
keep, buf, cur_len = [], [], 0
with open(div_fna) as fh:
    for line in fh:
        line = line.rstrip()
        if line.startswith(">"):
            if buf and cur_len >= 100:
                keep.extend(buf)
            buf = [line]
            cur_len = 0
        else:
            buf.append(line)
            cur_len += len(line)
    if buf and cur_len >= 100:
        keep.extend(buf)

with open(div_fil, "w") as fh:
    fh.write("\n".join(keep) + "\n")

n_markers = sum(1 for l in keep if l.startswith(">"))
print(f"  ✓ {n_markers} markers ≥100 bp kept for BLAST")

# ════════════════════════════════════════════════════════════════════════════
# STEP 2 — Build multi-strain BLAST database
# ════════════════════════════════════════════════════════════════════════════
print("\n[2] Building BLAST database from 23-strain panel")
combined_fna = f"{WORK_DIR}/panel_combined.fna"
blast_db     = f"{WORK_DIR}/panel_db"

with open(combined_fna, "w") as out:
    for sid, (label, ptype, fpath) in ALL_STRAINS.items():
        print(f"  adding {sid} ← {os.path.basename(fpath)}")
        with open(fpath) as fh:
            for line in fh:
                if line.startswith(">"):
                    # Prefix every contig with strain ID so subject is parseable
                    contig = line[1:].split()[0]
                    out.write(f">{sid}|{contig}\n")
                else:
                    out.write(line)

run(f"makeblastdb -in {combined_fna} -dbtype nucl -out {blast_db} -parse_seqids",
    "makeblastdb")

# ════════════════════════════════════════════════════════════════════════════
# STEP 3 — BLASTn markers vs panel
# ════════════════════════════════════════════════════════════════════════════
print("\n[3] BLASTn: divergent markers vs 23-strain panel")
blast_out = f"{WORK_DIR}/blast_panel.tsv"
run(
    f"blastn -query {div_fil} -db {blast_db} "
    f"-outfmt '6 qseqid sseqid pident length qlen qstart qend sstart send evalue bitscore' "
    f"-perc_identity 85 -qcov_hsp_perc 50 "
    f"-max_target_seqs 500 -num_threads 4 "
    f"-out {blast_out}",
    "blastn"
)

with open(blast_out) as fh:
    nhits = sum(1 for _ in fh)
print(f"  ✓ {nhits} BLAST hits")

# ════════════════════════════════════════════════════════════════════════════
# STEP 4 — Presence / absence matrix
# ════════════════════════════════════════════════════════════════════════════
print("\n[4] Computing presence/absence matrix")

cols = ["qseqid","sseqid","pident","length","qlen","qstart","qend","sstart","send","evalue","bitscore"]
blast_df = pd.read_csv(blast_out, sep="\t", header=None, names=cols)

# Parse strain ID from sseqid (format: STRAIN_ID|contig)
blast_df["strain"] = blast_df["sseqid"].str.split("|").str[0]

# Coverage filter: hit covers ≥50% of query
blast_df["qcov"] = (blast_df["qend"] - blast_df["qstart"] + 1) / blast_df["qlen"]
blast_df = blast_df[blast_df["qcov"] >= 0.5]

# All markers
all_markers = [l[1:].split("::")[0] if "::" in l else l[1:].strip()
               for l in keep if l.startswith(">")]

# Build presence set per strain
present = {sid: set() for sid in STRAIN_ORDER}
for _, row in blast_df.iterrows():
    sid = row["strain"]
    if sid in present:
        present[sid].add(row["qseqid"])

# Build matrix: 1 = present in strain, 0 = absent
PA = pd.DataFrame(0, index=all_markers, columns=STRAIN_ORDER, dtype=int)
for sid in STRAIN_ORDER:
    for mk in present[sid]:
        if mk in PA.index:
            PA.loc[mk, sid] = 1

print(f"  Marker presence (hit ≥85%id, ≥50% qcov) per strain:")
for sid in STRAIN_ORDER:
    pct = PA[sid].sum() / len(PA) * 100
    print(f"    {sid:<20} {PA[sid].sum():4d}/{len(PA)} markers present ({pct:.1f}%)")

# Sort markers: most pathogen-specific first (highest absence in commensals+outgroups)
non_path_cols  = list(COMMENSALS) + list(OUTGROUPS)
path_cols      = list(PATHOGENS)
PA["nonpath_hits"] = PA[non_path_cols].sum(axis=1)
PA["path_hits"]    = PA[path_cols].sum(axis=1)
PA_sorted = PA.sort_values(["nonpath_hits", "path_hits"], ascending=[True, False])
PA_sorted = PA_sorted.drop(columns=["nonpath_hits","path_hits"])

heatmap_data = PA_sorted[STRAIN_ORDER].values.astype(float)
n_markers_plot = len(PA_sorted)

# ════════════════════════════════════════════════════════════════════════════
# STEP 5 — Heatmap with row dendrogram
# ════════════════════════════════════════════════════════════════════════════
print("\n[5] Rendering heatmap with row dendrogram")

from scipy.spatial.distance import pdist
from scipy.cluster.hierarchy import linkage, dendrogram as sch_dendrogram  # dendrogram kept for leaf_order

# Colour: 0 = absent = cool steel-blue; 1 = present = warm red
cmap = ListedColormap(["#92C5DE", "#D6604D"])   # 0=absent→blue, 1=present→red

GROUP_COLOUR = {
    "EHEC":       "#1A7A6E",
    "UPEC":       "#E67E22",
    "ETEC":       "#8E44AD",
    "EPEC":       "#F39C12",
    "APEC":       "#D35400",
    "COMMENSAL":  "#2980B9",
    "K12":        "#1A5276",
    "LAB_B":      "#117A65",
    "PROBIOTIC":  "#1E8449",
    "OTHER_ESCH": "#7D3C98",
    "OTHER_GENUS":"#7F8C8D",
}

strain_labels = [ALL_STRAINS[s][0] for s in STRAIN_ORDER]
strain_ptypes = [ALL_STRAINS[s][1] for s in STRAIN_ORDER]

# ── row clustering: Jaccard distance + average linkage ───────────────────────
row_dist    = pdist(heatmap_data, metric="jaccard")
row_linkage = linkage(row_dist, method="average")

# get leaf order from dendrogram
dend_info   = sch_dendrogram(row_linkage, no_plot=True)
leaf_order  = dend_info["leaves"]
heatmap_data_ordered = heatmap_data[leaf_order, :]

# Summary bar: % markers present per strain (uses original sort, same values)
pct_present = [PA_sorted[s].sum() / n_markers_plot * 100 for s in STRAIN_ORDER]

# ── layout: 2-row × 1-col GridSpec (no dendrogram panel) ────────────────────
fig = plt.figure(figsize=(16, 11))
gs  = gridspec.GridSpec(
    2, 1,
    height_ratios=[1, 10],
    hspace=0.06,
    left=0.06, right=0.92,
    top=0.91, bottom=0.20
)

# ── bar chart ─────────────────────────────────────────────────────────────────
ax_bar = fig.add_subplot(gs[0, 0])
bar_colours = [GROUP_COLOUR[pt] for pt in strain_ptypes]
ax_bar.bar(range(len(STRAIN_ORDER)), pct_present,
           color=bar_colours, edgecolor="white", linewidth=0.5)
ax_bar.set_xlim(-0.5, len(STRAIN_ORDER)-0.5)
ax_bar.set_ylim(0, 140)
ax_bar.set_xticks([])
ax_bar.set_ylabel("% markers\npresent", fontsize=7, labelpad=4)
ax_bar.yaxis.set_major_locator(MaxNLocator(nbins=3, integer=True))
ax_bar.tick_params(axis="y", labelsize=6)
ax_bar.spines[["top","right","bottom"]].set_visible(False)
ax_bar.axhline(50, color="grey", linewidth=0.6, linestyle="--", alpha=0.5)
for xmid, label in [(4.5,"Pathogenic E. coli"),(14.5,"Commensal / Lab E. coli"),(21.0,"Outgroups")]:
    ax_bar.text(xmid, 135, label, ha="center", va="top", fontsize=7.5,
                fontweight="bold", color="black")

# ── main heatmap ──────────────────────────────────────────────────────────────
ax_main = fig.add_subplot(gs[1, 0])
ax_main.imshow(
    heatmap_data_ordered, aspect="auto", cmap=cmap,
    vmin=0, vmax=1, interpolation="none"
)

ax_main.set_xticks(range(len(STRAIN_ORDER)))
ax_main.set_xticklabels(strain_labels, rotation=45, ha="right", fontsize=7.5)
for tick, pt in zip(ax_main.get_xticklabels(), strain_ptypes):
    tick.set_color(GROUP_COLOUR.get(pt, "black"))

ax_main.set_yticks([])
ax_main.set_ylabel(f"Divergent markers (n={n_markers_plot})", fontsize=8, labelpad=4)

# vertical group dividers
for x in [9.5, 19.5]:
    ax_main.axvline(x=x, color="white", linewidth=1.5, linestyle="--")
    ax_bar.axvline( x=x, color="white", linewidth=1.5, linestyle="--")

# ── title ─────────────────────────────────────────────────────────────────────
fig.suptitle(
    "Distribution of Sakai O157:H7 Pilot Markers across E. coli and Outgroup Panel\n"
    f"(Sakai vs SE11 divergent regions, n={n_markers_plot} markers ≥100 bp | "
    "BLASTn ≥85% id, ≥50% qcov · rows clustered by Jaccard / average linkage)",
    fontsize=10, fontweight="bold", y=0.98
)

# ── legends ───────────────────────────────────────────────────────────────────
signal_entries = [
    mpatches.Patch(color="#D6604D", label="Present  (BLASTn hit ≥85% id, ≥50% qcov)"),
    mpatches.Patch(color="#92C5DE", label="Absent   (Sakai-specific, not found in strain)"),
]
pathotype_entries = [
    mpatches.Patch(color="#1A7A6E", label="EHEC (O157:H7 and non-O157)"),
    mpatches.Patch(color="#E67E22", label="UPEC"),
    mpatches.Patch(color="#8E44AD", label="ETEC"),
    mpatches.Patch(color="#F39C12", label="EPEC"),
    mpatches.Patch(color="#D35400", label="APEC"),
    mpatches.Patch(color="#2980B9", label="Commensal E. coli"),
    mpatches.Patch(color="#1A5276", label="K-12 laboratory strains"),
    mpatches.Patch(color="#117A65", label="Lab-B / Probiotic"),
    mpatches.Patch(color="#7D3C98", label="Other Escherichia spp."),
    mpatches.Patch(color="#7F8C8D", label="Other genus (S. aureus)"),
]

leg1 = fig.legend(handles=signal_entries,
                  title="Signal", title_fontsize=7,
                  loc="lower center", bbox_to_anchor=(0.28, 0.06),
                  ncol=2, fontsize=6.5, frameon=True, framealpha=0.8,
                  handlelength=1.2)
fig.legend(handles=pathotype_entries,
           title="Strain category", title_fontsize=7,
           loc="lower center", bbox_to_anchor=(0.65, 0.06),
           ncol=5, fontsize=6.5, frameon=True, framealpha=0.8,
           handlelength=1.2)
fig.add_artist(leg1)

out_path = f"{FIG_DIR}/fig_pilot_marker_heatmap.png"
plt.savefig(out_path, dpi=200, bbox_inches="tight", facecolor="white")
print(f"\n  ✓ Saved → {out_path}")
plt.close()

# ── save presence/absence table ───────────────────────────────────────────────
pa_out = f"{WORK_DIR}/pilot_presence_absence.tsv"
PA_sorted[STRAIN_ORDER].to_csv(pa_out, sep="\t")
print(f"  ✓ PA matrix → {pa_out}")
print("\nDone.")
