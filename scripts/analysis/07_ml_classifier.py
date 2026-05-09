#!/usr/bin/env python3
"""
Phase 2 (Proof-of-Concept) — XGBoost classifier with 10-feature vector.

Feature vector per candidate genomic region:
  1. blastn_identity        — top BLASTn % identity against metagenome reads
  2. cai_score              — Codon Adaptation Index (Sharp & Li 1987; O157:H7 Sakai ORFeome)
  3. gc_delta               — |GC% of marker − mean GC% of O157:H7 core genome| (50.5%)
  4. srna_density           — AU-rich 10-mer fraction (sRNA binding site proxy; Peer & Margalit 2011)
  5. align_coverage         — fraction of marker length non-redundantly covered by BLAST hits
  6. kmer_deviation         — cosine distance of marker 4-mer profile from Sakai genome profile
  7. presence_pathogenic    — fraction of pathogenic strains covering the marker (NUCmer)
  8. presence_non_pathogenic — fraction of non-pathogenic strains covering the marker (NUCmer)
  9. pangenome_score        — presence_pathogenic − presence_non_pathogenic
  10. anvio_cluster_score    — mean Anvi'o gene-cluster enrichment score for overlapping Sakai genes

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
import os
import re
import subprocess
import sys
import warnings
from math import exp, log as mlog
from pathlib import Path
from itertools import product

import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio.Data import CodonTable
from Bio.SeqUtils import gc_fraction
import xgboost as xgb
import shap
from sklearn.model_selection import StratifiedGroupKFold, cross_validate
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

BIN_DIR = Path(sys.executable).parent

REPO_ROOT   = Path(__file__).parents[2]
MARKERS_DIR    = REPO_ROOT / "data" / "results" / "markers"
BLAST_DIR      = REPO_ROOT / "data" / "results" / "blast_screen"
ML_DIR         = REPO_ROOT / "data" / "results" / "ml"
MINIMAP2_DIR   = REPO_ROOT / "data" / "results" / "minimap2"
ML_DIR.mkdir(parents=True, exist_ok=True)
VIZ_DIR     = REPO_ROOT / "data" / "results" / "figures"
VIZ_DIR.mkdir(parents=True, exist_ok=True)

REFERENCE        = REPO_ROOT / "data" / "genomes" / "O157H7_Sakai" / "O157H7_Sakai.fna"
REFERENCE_REFORM = REPO_ROOT / "data" / "results" / "pangenome" / "reformatted" / "O157H7_Sakai.fna"
GENE_CLUSTERS    = REPO_ROOT / "data" / "results" / "pangenome" / "gene_clusters.tsv"
GENOME_MANIFEST  = REPO_ROOT / "data" / "genomes" / "genome_manifest.tsv"
PANGENOME_DIR    = REPO_ROOT / "data" / "results" / "pangenome"

TIER_ENCODING = {"CONSERVED": 0, "MODERATE": 1, "DIVERGED": 2}

# O157:H7 Sakai core genome mean GC%: 50.5% (Hayashi et al. 2001, PMID 11206551)
REFERENCE_GC = 50.5


# ── Feature computation ───────────────────────────────────────────────────────

def compute_gc_delta(seq: str) -> float:
    """Absolute difference between marker GC% and reference core GC%."""
    gc = gc_fraction(seq) * 100
    return abs(gc - REFERENCE_GC)


def _predict_sakai_cds(genome_fasta: Path, out_dir: Path) -> Path:
    """Run prodigal on Sakai genome to predict CDS; cache result to out_dir."""
    cds_fasta = out_dir / "sakai_cds.fna"
    if cds_fasta.exists():
        return cds_fasta
    prodigal = BIN_DIR / "prodigal"
    if not prodigal.exists():
        raise FileNotFoundError(
            "prodigal not found in fp_pipeline env — run: conda install -c bioconda prodigal"
        )
    log.info("Predicting Sakai CDS with prodigal (one-time setup) → %s", cds_fasta.name)
    subprocess.run(
        [str(prodigal), "-i", str(genome_fasta), "-d", str(cds_fasta),
         "-p", "single", "-q"],
        check=True, capture_output=True,
    )
    return cds_fasta


def _build_rscu(cds_fasta: Path) -> dict[str, float]:
    """
    Compute Relative Synonymous Codon Usage (RSCU) from CDS sequences.
    RSCU_ij = X_ij / (1/n_i * sum_j X_ij)  where n_i = synonymous family size.
    Codons in a synonymous family with zero total usage receive RSCU=1 (neutral).
    """
    table = CodonTable.unambiguous_dna_by_name["Standard"]
    aa_to_codons: dict[str, list[str]] = {}
    for codon, aa in table.forward_table.items():
        aa_to_codons.setdefault(aa, []).append(codon)

    codon_counts: dict[str, int] = {}
    for rec in SeqIO.parse(str(cds_fasta), "fasta"):
        seq = str(rec.seq).upper()
        for i in range(0, len(seq) - 2, 3):
            codon = seq[i:i + 3]
            if len(codon) == 3 and all(c in "ACGT" for c in codon):
                codon_counts[codon] = codon_counts.get(codon, 0) + 1

    rscu: dict[str, float] = {}
    for aa, codons in aa_to_codons.items():
        total = sum(codon_counts.get(c, 0) for c in codons)
        n = len(codons)
        for codon in codons:
            rscu[codon] = (codon_counts.get(codon, 0) * n / total) if total > 0 else 1.0
    return rscu


def _compute_w(rscu: dict[str, float]) -> dict[str, float]:
    """Relative adaptedness w_i = RSCU_i / max(RSCU_j) for synonymous family j."""
    table = CodonTable.unambiguous_dna_by_name["Standard"]
    aa_to_codons: dict[str, list[str]] = {}
    for codon, aa in table.forward_table.items():
        aa_to_codons.setdefault(aa, []).append(codon)

    w: dict[str, float] = {}
    for aa, codons in aa_to_codons.items():
        max_rscu = max((rscu.get(c, 0.0) for c in codons), default=1.0) or 1.0
        for codon in codons:
            w[codon] = rscu.get(codon, 0.0) / max_rscu
    return w


_CAI_CACHE: dict[str, dict] = {}


def compute_cai(seq: str) -> float:
    """
    Codon Adaptation Index (Sharp & Li 1987) relative to the O157:H7 Sakai ORFeome
    predicted by prodigal. Returns geometric mean of per-codon relative adaptedness
    (w_i) values. Returns 0.5 for sequences with fewer than 5 usable sense codons.
    """
    if "w" not in _CAI_CACHE:
        cds_fasta = _predict_sakai_cds(REFERENCE, ML_DIR)
        rscu = _build_rscu(cds_fasta)
        _CAI_CACHE["w"] = _compute_w(rscu)
        _CAI_CACHE["stops"] = set(CodonTable.unambiguous_dna_by_name["Standard"].stop_codons)
        log.info("CAI reference built from Sakai ORFeome: %d codons in w-table",
                 len(_CAI_CACHE["w"]))

    w = _CAI_CACHE["w"]
    stops = _CAI_CACHE["stops"]

    seq = seq.upper()
    log_sum, n = 0.0, 0
    for i in range(0, len(seq) - 2, 3):
        codon = seq[i:i + 3]
        if len(codon) == 3 and codon not in stops and codon in w:
            wi = w[codon]
            if wi > 0:
                log_sum += mlog(wi)
                n += 1
    return exp(log_sum / n) if n >= 5 else 0.5


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


def _merge_intervals(starts: list[int], ends: list[int]) -> int:
    """Total non-overlapping query length covered by (start, end) intervals."""
    if not starts:
        return 0
    intervals = sorted(zip(starts, ends))
    lo, hi = intervals[0]
    total = 0
    for s, e in intervals[1:]:
        if s <= hi:
            hi = max(hi, e)
        else:
            total += hi - lo + 1
            lo, hi = s, e
    total += hi - lo + 1
    return total


def compute_blast_features(marker_id: str, blast_df: pd.DataFrame) -> tuple[float, float]:
    """
    For a given marker, return (max_identity, coverage_fraction) from BLAST hits.
    Coverage uses interval merging on query coordinates to avoid double-counting
    overlapping hits from the same genomic region.
    """
    hits = blast_df[blast_df["qseqid"].str.startswith(marker_id)]
    if hits.empty:
        return 0.0, 0.0
    max_identity = hits["pident"].max()
    match = re.search(r"len=(\d+)", marker_id)
    marker_len = int(match.group(1)) if match else 1000
    covered = _merge_intervals(hits["qstart"].tolist(), hits["qend"].tolist())
    coverage = min(covered / max(marker_len, 1), 1.0)
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


def compute_anvio_cluster_score(meta: pd.DataFrame) -> pd.DataFrame:
    """
    Pathogen enrichment score derived from Anvi'o gene cluster presence/absence.

    Steps:
      1. Build reformatted→original contig name map for O157H7_Sakai (by FASTA order).
      2. Export Sakai gene calls from the contigs DB (or reuse cached TSV).
      3. For each gene cluster, compute fraction of pathogenic vs non-pathogenic
         genomes that carry it.
      4. For each marker, average the cluster enrichment scores of all Sakai genes
         overlapping the marker region.

    Returns DataFrame indexed by marker_id with column 'anvio_cluster_score'.
    """
    PATHOGENIC = {"EHEC", "UPEC", "ETEC", "EAEC", "EPEC", "NMEC", "AIEC"}

    # ── 1. Contig name mapping: reformatted → original ────────────────────────
    orig_ids     = [r.id for r in SeqIO.parse(str(REFERENCE), "fasta")]
    reformed_ids = [r.id for r in SeqIO.parse(str(REFERENCE_REFORM), "fasta")]
    reformed_to_orig = dict(zip(reformed_ids, orig_ids))

    # ── 2. Sakai gene calls with positions ────────────────────────────────────
    sakai_db     = PANGENOME_DIR / "contigs_dbs" / "O157H7_Sakai.db"
    gene_calls_f = PANGENOME_DIR / "sakai_gene_calls.tsv"
    if not gene_calls_f.exists():
        anvi_bin = Path(
            subprocess.check_output(
                ["which", "anvi-export-gene-calls"],
                env={**__import__("os").environ,
                     "PATH": "/opt/anaconda3/envs/anvio8/bin:" + __import__("os").environ["PATH"]},
            ).decode().strip()
        )
        subprocess.run(
            [str(anvi_bin), "-c", str(sakai_db), "--gene-caller", "prodigal",
             "-o", str(gene_calls_f)],
            env={**__import__("os").environ,
                 "PATH": "/opt/anaconda3/envs/anvio8/bin:" + __import__("os").environ["PATH"]},
            check=True, capture_output=True,
        )
    gene_calls = pd.read_csv(gene_calls_f, sep="\t",
                             usecols=["gene_callers_id", "contig", "start", "stop"])
    gene_calls["orig_contig"] = gene_calls["contig"].map(reformed_to_orig)

    # ── 3. Gene cluster membership ────────────────────────────────────────────
    gc = pd.read_csv(GENE_CLUSTERS, sep="\t",
                     usecols=["gene_caller_id", "gene_cluster_id", "genome_name"])
    sakai_gc = gc[gc["genome_name"] == "O157H7_Sakai"][["gene_caller_id", "gene_cluster_id"]]
    gene_calls = gene_calls.merge(sakai_gc, left_on="gene_callers_id",
                                  right_on="gene_caller_id", how="left")

    # ── 4. Pathogen enrichment per cluster ───────────────────────────────────
    manifest = pd.read_csv(GENOME_MANIFEST, sep="\t")
    pathogen_genomes = set(manifest.loc[manifest["pathotype"].isin(PATHOGENIC), "label"])
    nonpath_genomes  = set(manifest.loc[~manifest["pathotype"].isin(PATHOGENIC), "label"])
    n_path    = max(len(pathogen_genomes), 1)
    n_nonpath = max(len(nonpath_genomes), 1)

    # gene_clusters.tsv uses G_-prefixed safe labels for digit-starting genomes;
    # normalise to match manifest labels before membership lookup
    def _norm(name: str) -> str:
        return name[2:] if name.startswith("G_") else name

    cluster_genomes = gc.groupby("gene_cluster_id")["genome_name"].apply(
        lambda s: {_norm(n) for n in s}
    ).to_dict()

    cluster_scores: dict[str, float] = {}
    for cid, genomes in cluster_genomes.items():
        pres_path    = len(genomes & pathogen_genomes)    / n_path
        pres_nonpath = len(genomes & nonpath_genomes)     / n_nonpath
        cluster_scores[cid] = pres_path - pres_nonpath

    # ── 5. Per-marker score ───────────────────────────────────────────────────
    rows = []
    for _, m in meta.iterrows():
        contig  = m["contig"]
        m_start = int(m["start"])
        m_end   = int(m["end"])

        overlap = gene_calls[
            (gene_calls["orig_contig"] == contig) &
            (gene_calls["start"]       <= m_end)  &
            (gene_calls["stop"]        >= m_start) &
            gene_calls["gene_cluster_id"].notna()
        ]

        if overlap.empty:
            score = 0.0
        else:
            score = float(np.mean(
                [cluster_scores.get(cid, 0.0) for cid in overlap["gene_cluster_id"]]
            ))

        rows.append({"marker_id": m["marker_id"], "anvio_cluster_score": round(score, 4)})

    return pd.DataFrame(rows).set_index("marker_id")


def compute_cog_enrichment_score(meta: pd.DataFrame) -> pd.DataFrame:
    """
    Per-marker COG functional enrichment score derived from
    anvi-compute-functional-enrichment-in-pan output.

    For each significantly PATHOGEN-enriched COG function (q < 0.05), the
    enrichment score is assigned to every gene cluster carrying that function.
    For each marker, the score is the maximum enrichment score across all
    overlapping Sakai gene clusters.  Normalised to [0, 1] by dividing by the
    observed max score across all clusters.
    """
    enrichment_file = PANGENOME_DIR / "enrichment_scores.tsv"
    if not enrichment_file.exists():
        log.warning("enrichment_scores.tsv not found — cog_enrichment_score will be 0")
        return pd.DataFrame(
            [{"marker_id": m, "cog_enrichment_score": 0.0} for m in meta["marker_id"]]
        ).set_index("marker_id")

    enrich = pd.read_csv(enrichment_file, sep="\t")
    sig = enrich[
        (enrich["adjusted_q_value"] < 0.05) &
        (enrich["associated_groups"] == "PATHOGEN")
    ]

    # Build cluster → enrichment score map (max score when cluster appears in
    # multiple significant functions — rare but possible)
    cluster_enrich: dict[str, float] = {}
    for _, row in sig.iterrows():
        score = float(row["enrichment_score"])
        for cid in str(row["gene_clusters_ids"]).split(","):
            cid = cid.strip()
            if cid:
                cluster_enrich[cid] = max(cluster_enrich.get(cid, 0.0), score)

    max_score = max(cluster_enrich.values()) if cluster_enrich else 1.0

    # Reuse Sakai gene→cluster mapping from the cached gene calls
    gene_calls_f = PANGENOME_DIR / "sakai_gene_calls.tsv"
    if not gene_calls_f.exists():
        log.warning("sakai_gene_calls.tsv not found — cog_enrichment_score will be 0")
        return pd.DataFrame(
            [{"marker_id": m, "cog_enrichment_score": 0.0} for m in meta["marker_id"]]
        ).set_index("marker_id")

    orig_ids     = [r.id for r in SeqIO.parse(str(REFERENCE), "fasta")]
    reformed_ids = [r.id for r in SeqIO.parse(str(REFERENCE_REFORM), "fasta")]
    reformed_to_orig = dict(zip(reformed_ids, orig_ids))

    gene_calls = pd.read_csv(gene_calls_f, sep="\t",
                             usecols=["gene_callers_id", "contig", "start", "stop"])
    gene_calls["orig_contig"] = gene_calls["contig"].map(reformed_to_orig)

    gc = pd.read_csv(GENE_CLUSTERS, sep="\t",
                     usecols=["gene_caller_id", "gene_cluster_id", "genome_name"])
    sakai_gc = gc[gc["genome_name"] == "O157H7_Sakai"][["gene_caller_id", "gene_cluster_id"]]
    gene_calls = gene_calls.merge(sakai_gc, left_on="gene_callers_id",
                                  right_on="gene_caller_id", how="left")

    rows = []
    for _, m in meta.iterrows():
        contig  = m["contig"]
        m_start = int(m["start"])
        m_end   = int(m["end"])

        overlap = gene_calls[
            (gene_calls["orig_contig"] == contig) &
            (gene_calls["start"]       <= m_end)  &
            (gene_calls["stop"]        >= m_start) &
            gene_calls["gene_cluster_id"].notna()
        ]

        if overlap.empty:
            score = 0.0
        else:
            raw = max(cluster_enrich.get(cid, 0.0)
                      for cid in overlap["gene_cluster_id"])
            score = round(raw / max_score, 4)

        rows.append({"marker_id": m["marker_id"], "cog_enrichment_score": score})

    log.info("  COG enrichment: %d significant PATHOGEN-enriched functions, "
             "%d clusters with scores", len(sig), len(cluster_enrich))
    return pd.DataFrame(rows).set_index("marker_id")


def load_minimap2_features(scheme: str = "A", window_size: int = 500) -> pd.DataFrame:
    """Load minimap2 divergence gradient features if available (minimap2-gradient branch).

    Matches candidates to NUCmer markers by contig + coordinate overlap.
    Returns empty DataFrame if minimap2 results don't exist (falls back to core features).
    """
    candidates_path = MINIMAP2_DIR / f"candidates_scheme{scheme}_{window_size}bp.tsv"
    if not candidates_path.exists():
        log.info("minimap2 candidates not found — running without gradient features")
        return pd.DataFrame()

    cands = pd.read_csv(candidates_path, sep="\t")
    log.info("Loaded %d minimap2 candidates (scheme %s, %dbp)", len(cands), scheme, window_size)
    return cands


def _overlap_minimap2(
    meta_row,
    mm_cands: pd.DataFrame,
) -> dict:
    """Find the best-overlapping minimap2 candidate for a NUCmer marker."""
    contig  = meta_row["contig"]
    m_start = int(meta_row["start"])
    m_end   = int(meta_row["end"])

    overlap = mm_cands[
        (mm_cands["contig"]        == contig) &
        (mm_cands["region_start"]  <  m_end) &
        (mm_cands["region_end"]    >  m_start)
    ]

    if overlap.empty:
        return {
            "mean_divergence":          0.0,
            "flank_conservation_2000bp": 0.5,
            "marker_score_2000bp":       0.0,
            "proportion_high_windows":   0.0,
            "proportion_mid_windows":    0.0,
        }

    # Pick candidate with maximum overlap
    best = overlap.iloc[
        ((overlap["region_end"].clip(upper=m_end) -
          overlap["region_start"].clip(lower=m_start)).argmax())
    ]
    return {
        "mean_divergence":           float(best.get("mean_divergence", 0.0)),
        "flank_conservation_2000bp": float(best.get("flank_conservation_2000bp", 0.5)),
        "marker_score_2000bp":       float(best.get("marker_score_2000bp", 0.0)),
        "proportion_high_windows":   float(best.get("proportion_high_windows", 0.0)),
        "proportion_mid_windows":    float(best.get("proportion_mid_windows", 0.0)),
    }


def build_feature_matrix(
    meta: pd.DataFrame,
    blast_df: pd.DataFrame,
    ref_seqs: dict[str, str],
    pan_presence: pd.DataFrame,
    anvio_scores: pd.DataFrame,
    cog_scores: pd.DataFrame,
) -> pd.DataFrame:
    mm_cands = load_minimap2_features(scheme="A", window_size=500)
    has_minimap2 = not mm_cands.empty
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

        if mid in pan_presence.index:
            pan_score    = pan_presence.loc[mid, "pangenome_score"]
            pres_path    = pan_presence.loc[mid, "presence_pathogenic"]
            pres_nonpath = pan_presence.loc[mid, "presence_non_pathogenic"]
        else:
            pan_score    = 0.0
            pres_path    = 0.0
            pres_nonpath = 0.0

        anvio_score = float(anvio_scores.loc[mid, "anvio_cluster_score"]) \
            if mid in anvio_scores.index else 0.0
        cog_score = float(cog_scores.loc[mid, "cog_enrichment_score"]) \
            if mid in cog_scores.index else 0.0

        mm_feats = _overlap_minimap2(m, mm_cands) if has_minimap2 else {
            "mean_divergence":           0.0,
            "flank_conservation_2000bp": 0.5,
            "marker_score_2000bp":       0.0,
            "proportion_high_windows":   0.0,
            "proportion_mid_windows":    0.0,
        }

        rows.append({
            "marker_id":               mid,
            "contig":                  contig,
            "start":                   start,
            "end":                     end,
            "tier":                    tier,
            "blastn_identity":         max_id,
            "cai_score":               compute_cai(seq),
            "gc_delta":                compute_gc_delta(seq),
            "srna_density":            compute_srna_density_proxy(seq),
            "align_coverage":          coverage,
            "kmer_deviation":          compute_kmer_deviation(seq, ref_seqs),
            "presence_pathogenic":     pres_path,
            "presence_non_pathogenic": pres_nonpath,
            "pangenome_score":         pan_score,
            "anvio_cluster_score":     anvio_score,
            "cog_enrichment_score":    cog_score,
            **mm_feats,
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
    "anvio_cluster_score",
    # minimap2 divergence gradient features (branch: minimap2-gradient)
    # populated when data/results/minimap2/candidates_schemeA_500bp.tsv exists
    "mean_divergence", "flank_conservation_2000bp", "marker_score_2000bp",
    "proportion_high_windows", "proportion_mid_windows",
]

# Features always present (NUCmer pipeline); minimap2 features are optional
_CORE_FEATURES = [
    "blastn_identity", "cai_score", "gc_delta",
    "srna_density", "align_coverage",
    "kmer_deviation",
    "presence_pathogenic", "presence_non_pathogenic", "pangenome_score",
    "anvio_cluster_score",
]


def train_and_evaluate(df: pd.DataFrame) -> xgb.XGBClassifier:
    # Use all FEATURE_COLS that are actually populated (minimap2 cols may be zeros)
    available = [c for c in FEATURE_COLS if c in df.columns and df[c].notna().all()]
    # Drop minimap2 features if all zero (pipeline ran without 02b)
    mm_cols = ["mean_divergence", "flank_conservation_2000bp", "marker_score_2000bp",
                "proportion_high_windows", "proportion_mid_windows"]
    has_mm = any(df.get(c, pd.Series([0])).sum() > 0 for c in mm_cols)
    if not has_mm:
        available = [c for c in available if c not in mm_cols]
    log.info("Training with %d features%s", len(available),
             " (+ minimap2 gradient features)" if has_mm else "")
    X = df[available].values
    y = df["label"].values
    # Group by genomic bin (500 kb windows on the main chromosome) so that
    # spatially adjacent markers are never split across train/test folds.
    # Using raw contig gives only 3 groups (2 plasmids + chromosome), which
    # is insufficient for 5-fold CV. 500 kb bins on NC_002695.2 yield ~11
    # groups; plasmid contigs are each treated as a single group.
    BIN_SIZE = 500_000
    def _make_bin(row):
        if row["contig"] == "NC_002695.2":
            return f"NC_002695.2_bin{int(row['start']) // BIN_SIZE:02d}"
        return row["contig"]
    groups = df.apply(_make_bin, axis=1).values
    n_groups = len(set(groups))
    n_splits = min(5, n_groups)
    if n_splits < 5:
        log.warning("Only %d spatial groups — reducing CV to %d folds", n_groups, n_splits)

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

    # StratifiedGroupKFold: preserves label ratio across folds AND keeps all
    # markers from the same contig in the same fold, eliminating leakage from
    # spatially correlated markers (adjacent windows share NUCmer block context).
    cv = StratifiedGroupKFold(n_splits=n_splits)
    cv_results = cross_validate(
        clf, X, y, cv=cv, groups=groups,
        scoring=["roc_auc", "average_precision", "f1"],
        return_train_score=True,
    )

    log.info("%d-fold grouped CV (500kb bins) — AUROC: %.3f ± %.3f | AUPRC: %.3f ± %.3f | F1: %.3f ± %.3f",
             n_splits,
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
    meta     = load_marker_metadata()
    blast_df = load_blast_results()
    ref_seqs = {rec.id: str(rec.seq) for rec in SeqIO.parse(str(REFERENCE), "fasta")}

    tiered_blocks_path = REPO_ROOT / "data" / "results" / "tiered_blocks.tsv"
    log.info("Computing pangenome presence/absence from NUCmer blocks")
    pan_presence = compute_pangenome_presence(meta, tiered_blocks_path)
    log.info("  NUCmer pangenome score range: %.3f – %.3f",
             pan_presence["pangenome_score"].min(),
             pan_presence["pangenome_score"].max())

    log.info("Computing Anvi'o gene-cluster enrichment scores")
    anvio_scores = compute_anvio_cluster_score(meta)
    log.info("  Anvi'o cluster score range: %.3f – %.3f",
             anvio_scores["anvio_cluster_score"].min(),
             anvio_scores["anvio_cluster_score"].max())

    log.info("Computing COG14 functional enrichment scores")
    cog_scores = compute_cog_enrichment_score(meta)

    log.info("Building feature matrix for %d markers", len(meta))
    df = build_feature_matrix(meta, blast_df, ref_seqs, pan_presence, anvio_scores, cog_scores)

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
