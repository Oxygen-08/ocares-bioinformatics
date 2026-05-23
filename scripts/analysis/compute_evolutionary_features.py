#!/usr/bin/env python3
"""
Phase 2: Dual-Feature Engineering — Delta-CAI & RSCU Cosine Distance
=====================================================================
Loads the core baseline reference produced by calculate_core_reference.py
and computes two evolutionary features for every one of the 6,189 accessory
gene clusters (with first-pass priority on the 66 pathogen-enriched COG
clusters):

  • Cosine_Distance  — cosine distance between the cluster's 59-dim RSCU
                        vector and the core RSCU reference vector.
                        Range [0, 1]; values near 0 indicate the cluster
                        uses the same codon style as the core genome.
  • Delta_CAI        — mean(CAI_cluster) − median(CAI_core).
                        Negative → cluster translated less efficiently than
                        core; positive → more efficiently (unusual for HGT
                        candidates, which typically show negative Delta-CAI
                        shortly after horizontal acquisition).

Core median CAI is estimated from one representative gene per core cluster
(Sakai genome used preferentially; first alphabetical genome as fallback).

Outputs:
  data/metrics/accessory_evolutionary_features.tsv
    columns: gene_cluster_id, Cosine_Distance, Delta_CAI,
             n_instances, n_genomes, pathogen_enriched

Run inside the fp_pipeline conda environment:
    conda run -n fp_pipeline python scripts/analysis/compute_evolutionary_features.py
"""

import json
import logging
import sqlite3
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
from Bio.Seq import Seq

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s  %(levelname)-7s  %(message)s",
    datefmt="%H:%M:%S",
)
log = logging.getLogger(__name__)

REPO_ROOT    = Path(__file__).parents[2]
PAN_DIR      = REPO_ROOT / "data" / "results" / "pangenome"
CONTIGS_DBS  = PAN_DIR / "contigs_dbs"
GC_TSV       = PAN_DIR / "gene_clusters.tsv"
ENRICH_TSV   = PAN_DIR / "enrichment_scores.tsv"
REFERENCE_JSON = REPO_ROOT / "data" / "metrics" / "core_evolutionary_reference.json"
OUT_FILE     = REPO_ROOT / "data" / "metrics" / "accessory_evolutionary_features.tsv"

STOP_CODONS = frozenset({"TAA", "TAG", "TGA"})


# ── Utilities ─────────────────────────────────────────────────────────────────

def reverse_complement(seq: str) -> str:
    return str(Seq(seq).reverse_complement())


def _db_name(genome_name: str) -> str:
    """Resolve Anvi'o G_-prefixed numeric genome names to their DB filename."""
    if genome_name.startswith("G_") and genome_name[2:].isdigit():
        return genome_name[2:]
    return genome_name


# ── Sequence extraction ───────────────────────────────────────────────────────

def _fetch_genome_seqs(
    db_path: Path,
    gene_caller_ids: list[int],
) -> dict[int, str]:
    """
    Batch-fetch CDS nucleotide sequences from an Anvi'o contigs database.
    Returns {gene_caller_id: nt_sequence} for successfully extracted genes.
    """
    if not gene_caller_ids:
        return {}

    result: dict[int, str] = {}
    conn = sqlite3.connect(str(db_path))
    try:
        placeholders = ",".join("?" * len(gene_caller_ids))
        gene_rows = conn.execute(
            f"SELECT gene_callers_id, contig, start, stop, direction, partial "
            f"FROM genes_in_contigs WHERE gene_callers_id IN ({placeholders})",
            gene_caller_ids,
        ).fetchall()

        needed_contigs = {row[1] for row in gene_rows if row[5] == 0}
        if not needed_contigs:
            return {}

        contig_seqs = {
            row[0]: row[1]
            for row in conn.execute(
                f"SELECT contig, sequence FROM contig_sequences "
                f"WHERE contig IN ({','.join('?' * len(needed_contigs))})",
                list(needed_contigs),
            ).fetchall()
        }

        for gcid, contig_name, start, stop, direction, partial in gene_rows:
            if partial or contig_name not in contig_seqs:
                continue
            nt = contig_seqs[contig_name][start:stop].upper()
            if direction == "r":
                nt = reverse_complement(nt)
            if len(nt) >= 6 and len(nt) % 3 == 0:
                result[gcid] = nt
    finally:
        conn.close()
    return result


# ── CAI computation ───────────────────────────────────────────────────────────

def compute_cai(nt_seq: str, weights: dict[str, float]) -> float | None:
    """
    Codon Adaptation Index (Sharp & Li 1987): geometric mean of per-codon
    relative adaptiveness weights, excluding stop codons and unknown codons.
    Uses log-space averaging for numerical stability.
    Returns None if no valid codons found.
    """
    log_w: list[float] = []
    for i in range(0, len(nt_seq) - 3, 3):  # -3 to skip terminal stop
        codon = nt_seq[i:i + 3].upper()
        if codon in STOP_CODONS or codon not in weights:
            continue
        w = weights[codon]
        if w > 0:
            log_w.append(np.log(w))
    if not log_w:
        return None
    return float(np.exp(np.mean(log_w)))


# ── RSCU computation ──────────────────────────────────────────────────────────

def codon_counts_from_cds(nt_seq: str) -> dict[str, int]:
    """Count codons in a CDS (excluding terminal stop)."""
    counts: dict[str, int] = defaultdict(int)
    for i in range(0, len(nt_seq) - 3, 3):
        codon = nt_seq[i:i + 3].upper()
        if len(codon) == 3 and codon.isalpha():
            counts[codon] += 1
    return dict(counts)


def compute_rscu_vector(
    pooled_counts: dict[str, int],
    syn_groups: dict[str, list[str]],
    codon_order: list[str],
) -> np.ndarray:
    """
    Compute RSCU for the 59 non-excluded synonymous codons.
    Returns a numpy array in the deterministic order of codon_order.
    """
    rscu: dict[str, float] = {}
    for aa, codons in syn_groups.items():
        n = len(codons)
        total = sum(pooled_counts.get(c, 0) for c in codons)
        if total == 0:
            for c in codons:
                rscu[c] = 1.0  # uniform fallback
        else:
            exp = total / n
            for c in codons:
                rscu[c] = pooled_counts.get(c, 0) / exp
    return np.array([rscu.get(c, 1.0) for c in codon_order])


def cosine_distance(v1: np.ndarray, v2: np.ndarray) -> float:
    """Cosine distance in [0, 1]. Returns 1.0 if either vector is zero."""
    n1, n2 = np.linalg.norm(v1), np.linalg.norm(v2)
    if n1 == 0 or n2 == 0:
        return 1.0
    return float(1.0 - np.dot(v1, v2) / (n1 * n2))


# ── Core median CAI ───────────────────────────────────────────────────────────

def estimate_core_median_cai(
    core_ids: set[str],
    gc: pd.DataFrame,
    cai_weights: dict[str, float],
    preferred_genome: str = "O157H7_Sakai",
) -> float:
    """
    Compute the median CAI of one representative gene per core cluster,
    using the Sakai genome where available, else the first alphabetical genome.
    """
    log.info("Estimating core median CAI from %d core clusters…", len(core_ids))

    core_gc = gc[gc["gene_cluster_id"].isin(core_ids)].copy()

    # Pick one representative gene per cluster
    representatives: list[tuple[str, int]] = []  # (genome_name, gene_caller_id)
    for cluster_id, group in core_gc.groupby("gene_cluster_id"):
        sakai_rows = group[group["genome_name"] == preferred_genome]
        if not sakai_rows.empty:
            rep = sakai_rows.iloc[0]
        else:
            rep = group.sort_values("genome_name").iloc[0]
        representatives.append((rep["genome_name"], int(rep["gene_caller_id"]), cluster_id))

    # Group by genome to minimise DB connections
    by_genome: dict[str, list[int]] = defaultdict(list)
    gcid_to_genome = {}
    for gname, gcid, cid in representatives:
        by_genome[gname].append(gcid)
        gcid_to_genome[gcid] = gname

    cai_scores: list[float] = []
    for genome_name, gene_ids in by_genome.items():
        db_path = CONTIGS_DBS / f"{_db_name(genome_name)}.db"
        if not db_path.exists():
            continue
        seqs = _fetch_genome_seqs(db_path, gene_ids)
        for gcid, nt in seqs.items():
            score = compute_cai(nt, cai_weights)
            if score is not None:
                cai_scores.append(score)

    median_cai = float(np.median(cai_scores))
    log.info(
        "Core CAI stats — n=%d  median=%.4f  mean=%.4f  std=%.4f  min=%.4f  max=%.4f",
        len(cai_scores), median_cai, np.mean(cai_scores),
        np.std(cai_scores), np.min(cai_scores), np.max(cai_scores),
    )
    return median_cai


# ── Main ─────────────────────────────────────────────────────────────────────

def main() -> None:
    log.info("=== Phase 2: Dual-Feature Engineering — fp_pipeline conda environment ===")

    # ── Load reference ────────────────────────────────────────────────────────
    log.info("Loading core reference from %s", REFERENCE_JSON)
    with open(REFERENCE_JSON) as fh:
        ref = json.load(fh)

    rscu_ref_vec   = np.array([ref["rscu_reference_vector"][c] for c in ref["rscu_codon_order"]])
    cai_weights    = ref["cai_weight_matrix"]
    codon_order    = ref["rscu_codon_order"]    # deterministic 59-codon list

    # Rebuild synonymous groups for RSCU computation
    _GENETIC_CODE: dict[str, str] = {
        "TTT": "Phe", "TTC": "Phe",
        "TTA": "Leu", "TTG": "Leu", "CTT": "Leu", "CTC": "Leu", "CTA": "Leu", "CTG": "Leu",
        "ATT": "Ile", "ATC": "Ile", "ATA": "Ile",
        "GTT": "Val", "GTC": "Val", "GTA": "Val", "GTG": "Val",
        "TCT": "Ser", "TCC": "Ser", "TCA": "Ser", "TCG": "Ser", "AGT": "Ser", "AGC": "Ser",
        "CCT": "Pro", "CCC": "Pro", "CCA": "Pro", "CCG": "Pro",
        "ACT": "Thr", "ACC": "Thr", "ACA": "Thr", "ACG": "Thr",
        "GCT": "Ala", "GCC": "Ala", "GCA": "Ala", "GCG": "Ala",
        "TAT": "Tyr", "TAC": "Tyr",
        "CAT": "His", "CAC": "His",
        "CAA": "Gln", "CAG": "Gln",
        "AAT": "Asn", "AAC": "Asn",
        "AAA": "Lys", "AAG": "Lys",
        "GAT": "Asp", "GAC": "Asp",
        "GAA": "Glu", "GAG": "Glu",
        "TGT": "Cys", "TGC": "Cys",
        "CGT": "Arg", "CGC": "Arg", "CGA": "Arg", "CGG": "Arg", "AGA": "Arg", "AGG": "Arg",
        "GGT": "Gly", "GGC": "Gly", "GGA": "Gly", "GGG": "Gly",
    }
    syn_groups: dict[str, list[str]] = defaultdict(list)
    for codon, aa in _GENETIC_CODE.items():
        if codon in set(codon_order):
            syn_groups[aa].append(codon)
    syn_groups = dict(syn_groups)

    # ── Load gene cluster assignments ─────────────────────────────────────────
    log.info("Loading gene cluster assignments from %s", GC_TSV)
    gc = pd.read_csv(GC_TSV, sep="\t")
    n_total_genomes = gc["genome_name"].nunique()
    cluster_counts  = gc.groupby("gene_cluster_id")["genome_name"].nunique()
    core_ids        = set(cluster_counts[cluster_counts == n_total_genomes].index)
    accessory_ids   = set(cluster_counts[cluster_counts < n_total_genomes].index)
    log.info("Core: %d clusters | Accessory: %d clusters", len(core_ids), len(accessory_ids))

    # ── Load pathogen-enriched cluster IDs ───────────────────────────────────
    log.info("Loading enrichment scores for pathogen-enriched cluster tagging…")
    enrich = pd.read_csv(ENRICH_TSV, sep="\t")
    sig_path = enrich[
        (enrich["adjusted_q_value"] < 0.05) &
        (enrich["associated_groups"].str.contains("PATHOGEN", na=False))
    ]
    pathogen_enriched_clusters: set[str] = set()
    for _, row in sig_path.iterrows():
        for cid in str(row["gene_clusters_ids"]).split(","):
            pathogen_enriched_clusters.add(cid.strip())
    log.info("Pathogen-enriched clusters (66 COG functions): %d cluster IDs",
             len(pathogen_enriched_clusters))

    # ── Core median CAI ───────────────────────────────────────────────────────
    core_median_cai = estimate_core_median_cai(core_ids, gc, cai_weights)

    # ── Per-cluster feature accumulation ─────────────────────────────────────
    # cluster_codons : cluster_id → aggregated codon counts (across all instances)
    # cluster_cai    : cluster_id → list of per-gene CAI scores
    cluster_codons: dict[str, dict[str, int]] = defaultdict(lambda: defaultdict(int))
    cluster_cai:    dict[str, list[float]]    = defaultdict(list)

    accessory_gc = gc[gc["gene_cluster_id"].isin(accessory_ids)].copy()
    genomes      = sorted(accessory_gc["genome_name"].unique())

    log.info("Processing %d genomes for %d accessory clusters…",
             len(genomes), len(accessory_ids))

    for i, genome_name in enumerate(genomes, 1):
        db_path = CONTIGS_DBS / f"{_db_name(genome_name)}.db"
        if not db_path.exists():
            log.warning("[%d/%d] Missing DB for %s — skipping", i, len(genomes), genome_name)
            continue

        genome_rows  = accessory_gc[accessory_gc["genome_name"] == genome_name]
        gene_ids     = genome_rows["gene_caller_id"].astype(int).tolist()
        cluster_of   = dict(zip(
            genome_rows["gene_caller_id"].astype(int),
            genome_rows["gene_cluster_id"],
        ))

        seqs = _fetch_genome_seqs(db_path, gene_ids)

        for gcid, nt in seqs.items():
            cid = cluster_of[gcid]
            # Codon accumulation for RSCU
            for codon, cnt in codon_counts_from_cds(nt).items():
                cluster_codons[cid][codon] += cnt
            # Per-gene CAI for Delta-CAI
            score = compute_cai(nt, cai_weights)
            if score is not None:
                cluster_cai[cid].append(score)

        if i % 10 == 0 or i == len(genomes):
            log.info("[%d/%d] genomes done — clusters with data: %d",
                     i, len(genomes), len(cluster_codons))

    log.info("Sequence accumulation complete. Computing per-cluster features…")

    # ── Prioritise pathogen-enriched clusters in output ordering ─────────────
    # Report pathogen-enriched first, then remaining accessory clusters
    priority_ids    = sorted(pathogen_enriched_clusters & accessory_ids)
    remaining_ids   = sorted(accessory_ids - pathogen_enriched_clusters)
    ordered_ids     = priority_ids + remaining_ids

    rows: list[dict] = []
    n_no_data = 0

    for cluster_id in ordered_ids:
        is_enriched = cluster_id in pathogen_enriched_clusters
        codons      = cluster_codons.get(cluster_id)
        cai_list    = cluster_cai.get(cluster_id, [])

        if not codons:
            # Cluster is accessory but no sequences were extracted (all partial/missing)
            rows.append({
                "gene_cluster_id":  cluster_id,
                "Cosine_Distance":  1.0,   # neutral flag: maximum distance
                "Delta_CAI":        0.0,   # neutral flag: no deviation measurable
                "n_instances":      0,
                "n_genomes":        int(cluster_counts.get(cluster_id, 0)),
                "pathogen_enriched": int(is_enriched),
            })
            n_no_data += 1
            continue

        rscu_vec  = compute_rscu_vector(codons, syn_groups, codon_order)
        cos_dist  = cosine_distance(rscu_vec, rscu_ref_vec)
        mean_cai  = float(np.mean(cai_list)) if cai_list else None
        delta_cai = (mean_cai - core_median_cai) if mean_cai is not None else 0.0

        rows.append({
            "gene_cluster_id":  cluster_id,
            "Cosine_Distance":  round(cos_dist, 6),
            "Delta_CAI":        round(delta_cai, 6),
            "n_instances":      len(cai_list),
            "n_genomes":        int(cluster_counts.get(cluster_id, 0)),
            "pathogen_enriched": int(is_enriched),
        })

    out_df = pd.DataFrame(rows)
    log.info("Feature computation complete — %d clusters total, %d had no extractable sequences",
             len(out_df), n_no_data)

    # ── Summary statistics ────────────────────────────────────────────────────
    path_mask = out_df["pathogen_enriched"] == 1
    other_mask = ~path_mask

    def _stats(series: pd.Series, label: str) -> None:
        log.info(
            "%-45s  n=%d  mean=%.4f  var=%.4f  min=%.4f  max=%.4f",
            label, len(series), series.mean(), series.var(),
            series.min(), series.max(),
        )

    log.info("=== Phase 4 preview: distribution summary ===")
    _stats(out_df.loc[path_mask,  "Cosine_Distance"], "Cosine_Distance [pathogen-enriched]")
    _stats(out_df.loc[other_mask, "Cosine_Distance"], "Cosine_Distance [other accessory]")
    _stats(out_df.loc[path_mask,  "Delta_CAI"],       "Delta_CAI       [pathogen-enriched]")
    _stats(out_df.loc[other_mask, "Delta_CAI"],       "Delta_CAI       [other accessory]")

    # ── Save ──────────────────────────────────────────────────────────────────
    out_df.to_csv(OUT_FILE, sep="\t", index=False)
    log.info("Saved → %s", OUT_FILE)
    log.info("=== Phase 2 complete ===")


if __name__ == "__main__":
    main()
