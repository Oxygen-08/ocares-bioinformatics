#!/usr/bin/env python3
"""
Phase 1b — minimap2 Divergence Gradient Pipeline.

Replaces NUCmer block-tiering (02_nucmer_tiered.py) with a sliding-window
divergence score derived from minimap2 PAF alignments. Each Sakai genome
window receives a divergence_score = 1 - (coverage_fraction * identity_fraction)
aggregated across all 30 commensal comparisons. Four threshold schemes are run
in parallel and compared.

Scientific basis
────────────────
Pathogenicity islands are often large genomic regions present in pathogenic
strains and absent from related non-pathogens, frequently acquired by HGT and
compositionally distinct. (Hacker & Kaper 2000 PMID 11018140)

Comparative genomic island methods identify query-genome regions absent from
related genomes. (Langille et al. 2008 BMC Bioinformatics DOI 10.1186/1471-2105-9-329)

Minimap2 is a peer-reviewed pairwise nucleotide aligner suitable for genome-
scale assembly comparison. (Li 2018 DOI 10.1093/bioinformatics/bty191)

Claim framing (per sentinel-PI guardrail):
  The biological basis of differential presence, sequence divergence, and
  genomic-island detection is established. This implementation tests whether a
  three-state nucleotide divergence gradient can serve as a structured and
  interpretable feature representation for pathogen-marker discovery.

Outputs (data/results/minimap2/)
────────────────────────────────
  paf/                         — raw minimap2 PAF files per commensal
  window_scores_500bp.tsv      — per-window divergence across all commensals
  window_scores_1000bp.tsv
  candidates_schemeA_500bp.tsv — candidate regions per scheme/window-size
  candidates_schemeA_500bp.fna
  ... (A/B/C/D × 500/1000bp)
  threshold_comparison.tsv     — summary table across all schemes
  negative_control/            — pathogen-vs-pathogen comparison outputs
  minimap2_pipeline.log
"""

import argparse
import csv
import logging
import math
import subprocess
import sys
import unittest
from collections import defaultdict
from pathlib import Path
from typing import Dict, Iterator, List, Optional, Tuple

import numpy as np
import pandas as pd
from Bio import SeqIO

# ── Paths ─────────────────────────────────────────────────────────────────────
REPO_ROOT   = Path(__file__).resolve().parents[2]
GENOME_DIR  = REPO_ROOT / "data" / "genomes"
MANIFEST    = GENOME_DIR / "genome_manifest.tsv"
REFERENCE   = GENOME_DIR / "O157H7_Sakai" / "O157H7_Sakai.fna"
RESULTS_DIR = REPO_ROOT / "data" / "results" / "minimap2"
PAF_DIR     = RESULTS_DIR / "paf"
NEG_DIR     = RESULTS_DIR / "negative_control"
RESULTS_DIR.mkdir(parents=True, exist_ok=True)
PAF_DIR.mkdir(parents=True, exist_ok=True)
NEG_DIR.mkdir(parents=True, exist_ok=True)

# ── Parameters ────────────────────────────────────────────────────────────────
WINDOW_SIZES  = [500, 1000]
FLANK_SIZES   = [2_000, 5_000]
MIN_CANDIDATE = 500   # bp — discard candidates shorter than this
MINIMAP_PRESET = "asm5"   # same-species assembly alignment

PATHOGENIC_PATHOTYPES = {
    "EHEC", "UPEC", "ETEC", "EAEC", "EPEC", "NMEC", "AIEC",
}

# Four threshold schemes from the sentinel document
THRESHOLD_SCHEMES: Dict[str, Dict] = {
    "A": {"low_hi": 0.20, "mid_hi": 0.60, "label": "biological_operational"},
    "B": {"low_hi": 0.15, "mid_hi": 0.50, "label": "stricter_high"},
    "C": {"low_hi": 0.25, "mid_hi": 0.70, "label": "conservative_high"},
    "D": {"low_hi": None, "mid_hi": None, "label": "data_driven_tertiles"},
}

# ── Logging ───────────────────────────────────────────────────────────────────
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        logging.FileHandler(RESULTS_DIR / "minimap2_pipeline.log"),
        logging.StreamHandler(sys.stdout),
    ],
)
log = logging.getLogger(__name__)


# ═══════════════════════════════════════════════════════════════════════════════
# MODULE 1 — Manifest loading
# ═══════════════════════════════════════════════════════════════════════════════

def load_manifest(pathotypes: Optional[set] = None) -> pd.DataFrame:
    """Load genome manifest, optionally filtering to specific pathotypes.

    Pass pathotypes=None to get all; pass set() to get commensals only
    (i.e. rows whose pathotype is NOT in PATHOGENIC_PATHOTYPES).
    """
    df = pd.read_csv(MANIFEST, sep="\t")
    df = df[df["status"] == "OK"].copy()
    df = df[df["label"] != "O157H7_Sakai"].copy()   # exclude self-comparison
    if pathotypes is not None:
        df = df[df["pathotype"].isin(pathotypes)].copy()
    return df.reset_index(drop=True)


def get_commensal_genomes() -> pd.DataFrame:
    df = pd.read_csv(MANIFEST, sep="\t")
    df = df[df["status"] == "OK"].copy()
    df = df[df["label"] != "O157H7_Sakai"].copy()
    return df[~df["pathotype"].isin(PATHOGENIC_PATHOTYPES)].reset_index(drop=True)


def get_pathogenic_genomes() -> pd.DataFrame:
    df = pd.read_csv(MANIFEST, sep="\t")
    df = df[df["status"] == "OK"].copy()
    df = df[df["label"] != "O157H7_Sakai"].copy()
    return df[df["pathotype"].isin(PATHOGENIC_PATHOTYPES)].reset_index(drop=True)


# ═══════════════════════════════════════════════════════════════════════════════
# MODULE 2 — minimap2 alignment
# ═══════════════════════════════════════════════════════════════════════════════

def run_minimap2(
    query_fasta: Path,
    target_fasta: Path,
    output_paf: Path,
    preset: str = MINIMAP_PRESET,
) -> Path:
    """Align query against target using minimap2 and write PAF output.

    minimap2 -x asm5 target.fasta query.fasta > out.paf
    Query = Sakai (pathogen), target = commensal.
    PAF col 1 = query (Sakai contig), col 6 = target (commensal contig).
    """
    if output_paf.exists() and output_paf.stat().st_size > 0:
        log.info("  PAF exists, skipping: %s", output_paf.name)
        return output_paf

    cmd = [
        "minimap2", "-x", preset,
        "--cs",               # cs tag for detailed alignment info (optional)
        "-t", "4",            # 4 threads
        str(target_fasta),
        str(query_fasta),
    ]
    log.info("  minimap2: %s vs %s", query_fasta.name, target_fasta.name)
    with open(output_paf, "w") as fh:
        result = subprocess.run(cmd, stdout=fh, stderr=subprocess.PIPE, text=True)
    if result.returncode != 0:
        log.error("minimap2 failed for %s: %s", output_paf.name, result.stderr[-500:])
        output_paf.unlink(missing_ok=True)
        raise RuntimeError(f"minimap2 failed: {output_paf.name}")
    log.info("    → %s (%.1f KB)", output_paf.name, output_paf.stat().st_size / 1024)
    return output_paf


# ═══════════════════════════════════════════════════════════════════════════════
# MODULE 3 — PAF parsing
# ═══════════════════════════════════════════════════════════════════════════════

def parse_paf(paf_path: Path) -> pd.DataFrame:
    """Parse minimap2 PAF into a DataFrame.

    Relevant columns (0-indexed):
      0  query_name       Sakai contig
      1  query_length
      2  query_start
      3  query_end
      9  n_residue_matches   number of matching bases
      10 alignment_block_len total alignment block length (including gaps)
      11 mapping_quality
    """
    rows = []
    with open(paf_path) as fh:
        for line in fh:
            if not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 12:
                continue
            try:
                rows.append({
                    "query_name":   f[0],
                    "query_length": int(f[1]),
                    "query_start":  int(f[2]),
                    "query_end":    int(f[3]),
                    "n_matches":    int(f[9]),
                    "block_len":    int(f[10]),
                    "mapq":         int(f[11]),
                })
            except (ValueError, IndexError):
                continue
    if not rows:
        return pd.DataFrame(columns=[
            "query_name", "query_length", "query_start",
            "query_end", "n_matches", "block_len", "mapq",
        ])
    return pd.DataFrame(rows)


# ═══════════════════════════════════════════════════════════════════════════════
# MODULE 4 — Window scoring
# ═══════════════════════════════════════════════════════════════════════════════

def _merge_covered(starts: List[int], ends: List[int]) -> int:
    """Total non-overlapping bases covered by a set of (start, end) intervals."""
    if not starts:
        return 0
    pairs = sorted(zip(starts, ends))
    lo, hi = pairs[0]
    total = 0
    for s, e in pairs[1:]:
        if s <= hi:
            hi = max(hi, e)
        else:
            total += hi - lo
            lo, hi = s, e
    return total + (hi - lo)


def score_windows(
    paf_df: pd.DataFrame,
    contig_lengths: Dict[str, int],
    window_size: int,
) -> pd.DataFrame:
    """Compute divergence_score for every window across all Sakai contigs.

    divergence_score = 1 - (coverage_fraction * identity_fraction)

    coverage_fraction  = non-overlapping aligned bases in window / window_size
    identity_fraction  = weighted-average (n_matches / block_len) for overlapping
                         alignments, weighted by their overlap with the window.

    Windows with zero coverage → divergence_score = 1.0 (fully absent/divergent).
    """
    rows = []
    for contig, clen in contig_lengths.items():
        contig_hits = paf_df[paf_df["query_name"] == contig]
        for w_start in range(0, clen, window_size):
            w_end = min(w_start + window_size, clen)
            actual_size = w_end - w_start

            # Find alignments overlapping this window
            if not contig_hits.empty:
                overlap = contig_hits[
                    (contig_hits["query_start"] < w_end) &
                    (contig_hits["query_end"]   > w_start)
                ]
            else:
                overlap = pd.DataFrame()

            if overlap.empty:
                rows.append({
                    "contig":      contig,
                    "w_start":     w_start,
                    "w_end":       w_end,
                    "window_size": actual_size,
                    "cov_frac":    0.0,
                    "id_frac":     0.0,
                    "div_score":   1.0,
                    "n_hits":      0,
                })
                continue

            # Clip intervals to window boundaries
            clipped_starts = overlap["query_start"].clip(lower=w_start).tolist()
            clipped_ends   = overlap["query_end"].clip(upper=w_end).tolist()
            covered = _merge_covered(clipped_starts, clipped_ends)
            cov_frac = min(covered / actual_size, 1.0)

            # Weighted-average identity (weight = overlap length in window)
            weights = [e - s for s, e in zip(clipped_starts, clipped_ends)]
            identities = (overlap["n_matches"] / overlap["block_len"].clip(lower=1)).tolist()
            total_weight = sum(weights)
            if total_weight > 0:
                id_frac = sum(w * i for w, i in zip(weights, identities)) / total_weight
            else:
                id_frac = 0.0

            div_score = 1.0 - (cov_frac * id_frac)
            rows.append({
                "contig":      contig,
                "w_start":     w_start,
                "w_end":       w_end,
                "window_size": actual_size,
                "cov_frac":    round(cov_frac, 6),
                "id_frac":     round(id_frac, 6),
                "div_score":   round(div_score, 6),
                "n_hits":      len(overlap),
            })
    return pd.DataFrame(rows)


def aggregate_scores(all_scores: List[pd.DataFrame]) -> pd.DataFrame:
    """Aggregate per-comparison window scores across all commensal comparisons.

    Returns one row per window with:
      mean_div_score     — mean divergence across all comparisons
      min_div_score      — most conserved comparison (strictest)
      max_div_score      — most divergent comparison
      std_div_score      — variability across comparisons
      n_comparisons      — number of comparisons contributing
    """
    combined = pd.concat(all_scores, ignore_index=True)
    agg = (
        combined
        .groupby(["contig", "w_start", "w_end", "window_size"])["div_score"]
        .agg(
            mean_div_score="mean",
            min_div_score="min",
            max_div_score="max",
            std_div_score="std",
            n_comparisons="count",
        )
        .reset_index()
    )
    agg["std_div_score"] = agg["std_div_score"].fillna(0.0)
    return agg


# ═══════════════════════════════════════════════════════════════════════════════
# MODULE 5 — Gradient classification
# ═══════════════════════════════════════════════════════════════════════════════

def classify_gradient(
    scores: pd.Series,
    scheme: str,
    scheme_params: Dict,
) -> pd.Series:
    """Classify divergence scores into LOW / MID / HIGH.

    Scheme D uses empirical tertiles of the score distribution.
    """
    if scheme == "D":
        t1, t2 = scores.quantile([1/3, 2/3]).values
        low_hi, mid_hi = t1, t2
        log.info("  Scheme D tertiles: LOW≤%.3f, MID≤%.3f, HIGH>%.3f", t1, t2, t2)
    else:
        low_hi = scheme_params["low_hi"]
        mid_hi = scheme_params["mid_hi"]

    def _grade(v: float) -> str:
        if v <= low_hi:
            return "LOW"
        elif v <= mid_hi:
            return "MID"
        return "HIGH"

    return scores.map(_grade)


# ═══════════════════════════════════════════════════════════════════════════════
# MODULE 6 — Candidate region extraction
# ═══════════════════════════════════════════════════════════════════════════════

def extract_candidates(
    window_df: pd.DataFrame,
    gradient_col: str,
    min_length: int = MIN_CANDIDATE,
    bridge_mid: bool = True,
) -> pd.DataFrame:
    """Merge adjacent HIGH windows into candidate divergent regions.

    bridge_mid: allow bridging across exactly 1 MID window flanked by HIGH on both sides.
    Candidates shorter than min_length are discarded.
    """
    candidates = []
    for contig, grp in window_df.groupby("contig"):
        grp = grp.sort_values("w_start").reset_index(drop=True)
        grades = grp[gradient_col].tolist()
        starts = grp["w_start"].tolist()
        ends   = grp["w_end"].tolist()
        scores = grp["mean_div_score"].tolist()

        # Expand HIGH runs, optionally bridging single MID windows
        if bridge_mid:
            bridged = list(grades)
            for i in range(1, len(grades) - 1):
                if (grades[i] == "MID"
                        and grades[i - 1] == "HIGH"
                        and grades[i + 1] == "HIGH"):
                    bridged[i] = "HIGH"
            grades = bridged

        # Merge contiguous HIGH stretches
        in_region = False
        reg_start = reg_end = 0
        reg_scores: List[float] = []
        window_grades: List[str] = []

        def _close_region(rs, re, rsc, wg):
            length = re - rs
            if length >= min_length:
                low_p  = wg.count("LOW")  / max(len(wg), 1)
                mid_p  = wg.count("MID")  / max(len(wg), 1)
                high_p = wg.count("HIGH") / max(len(wg), 1)
                candidates.append({
                    "contig":                contig,
                    "region_start":          rs,
                    "region_end":            re,
                    "region_length":         length,
                    "mean_divergence":       round(float(np.mean(rsc)), 6),
                    "max_divergence":        round(float(np.max(rsc)), 6),
                    "median_divergence":     round(float(np.median(rsc)), 6),
                    "proportion_low_windows":  round(low_p, 4),
                    "proportion_mid_windows":  round(mid_p, 4),
                    "proportion_high_windows": round(high_p, 4),
                    "n_windows":             len(wg),
                })

        for i, grade in enumerate(grades):
            if grade == "HIGH":
                if not in_region:
                    in_region = True
                    reg_start = starts[i]
                    reg_scores = []
                    window_grades = []
                reg_end = ends[i]
                reg_scores.append(scores[i])
                window_grades.append(grade)
            else:
                if in_region:
                    _close_region(reg_start, reg_end, reg_scores, window_grades)
                    in_region = False
                    reg_scores = []
                    window_grades = []

        if in_region:
            _close_region(reg_start, reg_end, reg_scores, window_grades)

    return pd.DataFrame(candidates)


# ═══════════════════════════════════════════════════════════════════════════════
# MODULE 7 — Flanking conservation
# ═══════════════════════════════════════════════════════════════════════════════

def compute_flanking_conservation(
    candidates: pd.DataFrame,
    window_df: pd.DataFrame,
    flank_size: int,
) -> pd.DataFrame:
    """Compute flanking conservation score for each candidate region.

    Flanking windows upstream and downstream within flank_size bp are collected.
    mean_flank_divergence = mean divergence of those flanking windows.
    flank_conservation    = 1 - mean_flank_divergence.

    Rationale: a strongly divergent region embedded in conserved flanking context
    is biologically more consistent with an HGT-inserted genomic island than a
    random alignment artefact or a contig end.

    marker_score = mean_divergence * flank_conservation * log1p(region_length)
    """
    if candidates.empty:
        return candidates.copy()

    results = []
    for _, cand in candidates.iterrows():
        contig = cand["contig"]
        rs, re = int(cand["region_start"]), int(cand["region_end"])

        upstream   = window_df[
            (window_df["contig"] == contig) &
            (window_df["w_end"]   <= rs) &
            (window_df["w_start"] >= rs - flank_size)
        ]["mean_div_score"]

        downstream = window_df[
            (window_df["contig"] == contig) &
            (window_df["w_start"] >= re) &
            (window_df["w_end"]   <= re + flank_size)
        ]["mean_div_score"]

        flank_scores = pd.concat([upstream, downstream])
        mean_flank   = float(flank_scores.mean()) if not flank_scores.empty else 0.5
        flank_cons   = round(1.0 - mean_flank, 6)
        m_score      = round(
            float(cand["mean_divergence"]) * flank_cons * math.log1p(int(cand["region_length"])),
            6,
        )
        row = cand.to_dict()
        row[f"flank_conservation_{flank_size}bp"] = flank_cons
        row[f"marker_score_{flank_size}bp"]        = m_score
        results.append(row)

    return pd.DataFrame(results)


# ═══════════════════════════════════════════════════════════════════════════════
# MODULE 8 — Threshold validation summary
# ═══════════════════════════════════════════════════════════════════════════════

def build_threshold_comparison(
    scheme_results: Dict[str, Dict],
    contig_lengths: Dict[str, int],
) -> pd.DataFrame:
    """Build comparison table across all threshold schemes and window sizes."""
    genome_size = sum(contig_lengths.values())
    rows = []
    for key, res in scheme_results.items():
        scheme, ws = key.split("_")
        cands = res["candidates"]
        windows = res["windows"]
        high_windows = (windows["gradient"] == "HIGH").sum()
        rows.append({
            "scheme":                  scheme,
            "scheme_label":            THRESHOLD_SCHEMES[scheme]["label"],
            "window_size_bp":          int(ws.replace("bp", "")),
            "total_windows":           len(windows),
            "high_windows":            int(high_windows),
            "pct_genome_high":         round(100 * high_windows * int(ws.replace("bp", "")) / max(genome_size, 1), 2),
            "n_candidate_regions":     len(cands),
            "median_candidate_len_bp": int(cands["region_length"].median()) if not cands.empty else 0,
            "median_marker_score":     round(float(cands["marker_score_2000bp"].median()), 4) if not cands.empty and "marker_score_2000bp" in cands.columns else 0.0,
            "total_candidate_bp":      int(cands["region_length"].sum()) if not cands.empty else 0,
        })
    return pd.DataFrame(rows).sort_values(["window_size_bp", "scheme"])


# ═══════════════════════════════════════════════════════════════════════════════
# MODULE 9 — FASTA output
# ═══════════════════════════════════════════════════════════════════════════════

def write_candidate_fasta(
    candidates: pd.DataFrame,
    ref_seqs: Dict[str, str],
    output_fna: Path,
    scheme: str,
    window_size: int,
) -> int:
    """Write candidate regions as a FASTA file. Returns number of sequences written."""
    written = 0
    with open(output_fna, "w") as fh:
        for i, row in candidates.iterrows():
            contig = row["contig"]
            if contig not in ref_seqs:
                continue
            seq = ref_seqs[contig][int(row["region_start"]):int(row["region_end"])]
            if not seq:
                continue
            marker_id = f"MM_{scheme}_{i+1:04d}"
            ms = row.get("marker_score_2000bp", 0.0)
            fh.write(
                f">{marker_id} contig={contig} "
                f"start={row['region_start']} end={row['region_end']} "
                f"len={row['region_length']} "
                f"mean_div={row['mean_divergence']:.4f} "
                f"marker_score={ms:.4f} "
                f"scheme={scheme} window={window_size}bp\n"
            )
            fh.write(seq + "\n")
            written += 1
    return written


# ═══════════════════════════════════════════════════════════════════════════════
# MODULE 10 — Negative control: pathogen vs pathogen
# ═══════════════════════════════════════════════════════════════════════════════

def run_negative_control(
    ref_seqs: Dict[str, str],
    contig_lengths: Dict[str, int],
    pathogenic_df: pd.DataFrame,
    window_size: int = 500,
    max_comparisons: int = 5,
) -> Dict:
    """Compare Sakai against other pathogenic strains as a negative control.

    Biological expectation: fewer HIGH-divergence windows than pathogen-vs-commensal,
    because pathogens share more PAI/virulence sequence with each other.
    """
    log.info("Negative control: Sakai vs pathogenic strains (%d comparisons)", min(len(pathogenic_df), max_comparisons))
    sample = pathogenic_df.head(max_comparisons)
    all_scores = []
    for _, row in sample.iterrows():
        paf_path = NEG_DIR / f"sakai_vs_{row['label']}.paf"
        try:
            run_minimap2(REFERENCE, Path(row["fasta_path"]), paf_path)
            paf_df = parse_paf(paf_path)
            scores = score_windows(paf_df, contig_lengths, window_size)
            scores["comparison"] = row["label"]
            all_scores.append(scores)
        except Exception as e:
            log.warning("  Negative control failed for %s: %s", row["label"], e)

    if not all_scores:
        return {"error": "all negative control comparisons failed"}

    agg = aggregate_scores(all_scores)
    high_frac = (agg["mean_div_score"] > 0.60).mean()
    log.info("  Negative control: %.1f%% HIGH windows (pathogen-vs-pathogen)", 100 * high_frac)
    agg.to_csv(NEG_DIR / f"neg_control_window_scores_{window_size}bp.tsv", sep="\t", index=False)
    return {
        "n_comparisons":      len(all_scores),
        "pct_high_windows":   round(100 * high_frac, 2),
        "mean_div_score_mean": round(float(agg["mean_div_score"].mean()), 4),
    }


# ═══════════════════════════════════════════════════════════════════════════════
# MODULE 11 — Unit tests
# ═══════════════════════════════════════════════════════════════════════════════

class TestWindowScoring(unittest.TestCase):

    def _make_paf(self, rows):
        return pd.DataFrame(rows, columns=[
            "query_name", "query_length", "query_start",
            "query_end", "n_matches", "block_len", "mapq",
        ])

    def test_fully_aligned_window(self):
        """Window fully covered with perfect identity → divergence_score ≈ 0."""
        paf = self._make_paf([
            ("chr1", 10000, 0, 500, 500, 500, 60),
        ])
        result = score_windows(paf, {"chr1": 600}, window_size=500)
        row = result[result["w_start"] == 0].iloc[0]
        self.assertAlmostEqual(row["div_score"], 0.0, places=3)

    def test_partially_aligned_window(self):
        """Window 50% covered with 80% identity → divergence_score ≈ 1 - 0.5*0.8 = 0.6."""
        paf = self._make_paf([
            ("chr1", 10000, 0, 250, 200, 250, 60),  # 50% of 500bp window
        ])
        result = score_windows(paf, {"chr1": 600}, window_size=500)
        row = result[result["w_start"] == 0].iloc[0]
        expected = 1.0 - (0.5 * (200 / 250))
        self.assertAlmostEqual(row["div_score"], expected, places=3)

    def test_unaligned_window(self):
        """Window with no alignments → divergence_score = 1.0."""
        paf = self._make_paf([
            ("chr1", 10000, 1000, 1500, 500, 500, 60),  # alignment outside window
        ])
        result = score_windows(paf, {"chr1": 600}, window_size=500)
        row = result[result["w_start"] == 0].iloc[0]
        self.assertAlmostEqual(row["div_score"], 1.0, places=6)

    def test_adjacent_high_window_merging(self):
        """Two adjacent HIGH windows merge into one candidate region."""
        window_df = pd.DataFrame([
            {"contig": "c1", "w_start": 0,   "w_end": 500,  "mean_div_score": 0.9, "window_size": 500},
            {"contig": "c1", "w_start": 500, "w_end": 1000, "mean_div_score": 0.85, "window_size": 500},
            {"contig": "c1", "w_start": 1000,"w_end": 1500, "mean_div_score": 0.1, "window_size": 500},
        ])
        window_df["gradient"] = classify_gradient(window_df["mean_div_score"], "A", THRESHOLD_SCHEMES["A"])
        cands = extract_candidates(window_df, "gradient", min_length=500, bridge_mid=False)
        self.assertEqual(len(cands), 1)
        self.assertEqual(cands.iloc[0]["region_start"], 0)
        self.assertEqual(cands.iloc[0]["region_end"], 1000)

    def test_flank_conservation(self):
        """Region flanked by conserved windows → high flank_conservation."""
        window_df = pd.DataFrame([
            {"contig": "c1", "w_start": 0,    "w_end": 500,  "mean_div_score": 0.05},
            {"contig": "c1", "w_start": 500,  "w_end": 1000, "mean_div_score": 0.9},
            {"contig": "c1", "w_start": 1000, "w_end": 1500, "mean_div_score": 0.05},
        ])
        candidates = pd.DataFrame([{
            "contig": "c1", "region_start": 500, "region_end": 1000,
            "region_length": 500, "mean_divergence": 0.9,
            "max_divergence": 0.9, "median_divergence": 0.9,
            "proportion_low_windows": 0.0, "proportion_mid_windows": 0.0,
            "proportion_high_windows": 1.0, "n_windows": 1,
        }])
        result = compute_flanking_conservation(candidates, window_df, flank_size=2000)
        self.assertGreater(result.iloc[0]["flank_conservation_2000bp"], 0.8)

    def test_merge_covered(self):
        """Interval merging handles overlapping intervals correctly."""
        self.assertEqual(_merge_covered([0, 100, 150], [100, 200, 250]), 250)
        self.assertEqual(_merge_covered([0, 50], [50, 100]), 100)
        self.assertEqual(_merge_covered([], []), 0)


# ═══════════════════════════════════════════════════════════════════════════════
# MAIN
# ═══════════════════════════════════════════════════════════════════════════════

def main(run_tests: bool = False, window_sizes: Optional[List[int]] = None) -> None:

    if run_tests:
        log.info("Running unit tests…")
        suite = unittest.TestLoader().loadTestsFromTestCase(TestWindowScoring)
        runner = unittest.TextTestRunner(verbosity=2)
        result = runner.run(suite)
        if not result.wasSuccessful():
            sys.exit(1)
        log.info("All unit tests passed.")
        return

    ws_list = window_sizes or WINDOW_SIZES

    # ── Load reference genome ─────────────────────────────────────────────────
    log.info("Loading reference: %s", REFERENCE)
    ref_seqs      = {r.id: str(r.seq) for r in SeqIO.parse(REFERENCE, "fasta")}
    contig_lengths = {k: len(v) for k, v in ref_seqs.items()}
    log.info("  Contigs: %s", {k: v for k, v in contig_lengths.items()})

    # ── Load genome manifests ─────────────────────────────────────────────────
    commensal_df  = get_commensal_genomes()
    pathogenic_df = get_pathogenic_genomes()
    log.info("Commensal genomes: %d | Pathogenic genomes: %d",
             len(commensal_df), len(pathogenic_df))

    # ── Run minimap2 for each commensal ───────────────────────────────────────
    log.info("Step 1: Aligning Sakai vs %d commensal genomes…", len(commensal_df))
    all_paf_dfs = []
    for _, row in commensal_df.iterrows():
        label     = row["label"]
        fasta     = Path(row["fasta_path"])
        paf_path  = PAF_DIR / f"sakai_vs_{label}.paf"
        try:
            run_minimap2(REFERENCE, fasta, paf_path)
            paf_df = parse_paf(paf_path)
            paf_df["comparison"] = label
            all_paf_dfs.append(paf_df)
        except Exception as e:
            log.warning("  Skipping %s: %s", label, e)

    if not all_paf_dfs:
        log.error("No PAF files parsed — check minimap2 installation and genome paths.")
        sys.exit(1)

    # ── Score windows and aggregate ───────────────────────────────────────────
    scheme_results: Dict[str, Dict] = {}

    for ws in ws_list:
        log.info("Step 2: Scoring windows (%d bp) across %d comparisons…", ws, len(all_paf_dfs))
        per_comparison_scores = []
        for paf_df in all_paf_dfs:
            scores = score_windows(paf_df, contig_lengths, ws)
            per_comparison_scores.append(scores)

        agg_df = aggregate_scores(per_comparison_scores)
        # Also compute fraction of comparisons where this window is HIGH (scheme A threshold)
        high_counts = []
        for paf_df in per_comparison_scores:
            high_counts.append((paf_df["div_score"] > 0.60).astype(int).rename("high"))
        agg_df["frac_comparisons_high"] = sum(hc.values for hc in high_counts) / max(len(high_counts), 1)

        agg_df.to_csv(RESULTS_DIR / f"window_scores_{ws}bp.tsv", sep="\t", index=False)
        log.info("  Window scores saved: %d windows", len(agg_df))

        # ── Apply all threshold schemes ───────────────────────────────────────
        for scheme, params in THRESHOLD_SCHEMES.items():
            log.info("  Scheme %s (%s)…", scheme, params["label"])
            window_df = agg_df.copy()
            window_df["gradient"] = classify_gradient(
                window_df["mean_div_score"], scheme, params
            )

            high_count = (window_df["gradient"] == "HIGH").sum()
            pct_high   = 100 * high_count * ws / max(sum(contig_lengths.values()), 1)
            log.info("    HIGH windows: %d (%.1f%% of genome)", high_count, pct_high)

            # Extract candidate regions
            cands = extract_candidates(window_df, "gradient", min_length=MIN_CANDIDATE)
            log.info("    Candidate regions: %d", len(cands))

            if not cands.empty:
                # Flanking conservation for both flank sizes
                for fs in FLANK_SIZES:
                    cands = compute_flanking_conservation(cands, window_df, fs)

                out_tsv = RESULTS_DIR / f"candidates_scheme{scheme}_{ws}bp.tsv"
                cands.to_csv(out_tsv, sep="\t", index=False)

                out_fna = RESULTS_DIR / f"candidates_scheme{scheme}_{ws}bp.fna"
                n_written = write_candidate_fasta(cands, ref_seqs, out_fna, scheme, ws)
                log.info("    FASTA written: %d sequences → %s", n_written, out_fna.name)

            key = f"{scheme}_{ws}bp"
            scheme_results[key] = {"windows": window_df, "candidates": cands}

    # ── Threshold comparison table ────────────────────────────────────────────
    log.info("Step 3: Building threshold comparison table…")
    comp_df = build_threshold_comparison(scheme_results, contig_lengths)
    comp_df.to_csv(RESULTS_DIR / "threshold_comparison.tsv", sep="\t", index=False)
    log.info("\n%s", comp_df.to_string(index=False))

    # ── Negative control ──────────────────────────────────────────────────────
    log.info("Step 4: Negative control (Sakai vs pathogenic strains)…")
    neg_ctrl = run_negative_control(ref_seqs, contig_lengths, pathogenic_df)
    log.info("  Result: %s", neg_ctrl)

    # Save negative control summary
    with open(NEG_DIR / "negative_control_summary.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=neg_ctrl.keys(), delimiter="\t")
        w.writeheader()
        w.writerow(neg_ctrl)

    log.info("Pipeline complete. Results: %s", RESULTS_DIR)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--test", action="store_true", help="Run unit tests and exit")
    parser.add_argument("--window-sizes", nargs="+", type=int,
                        default=WINDOW_SIZES, metavar="N",
                        help="Window sizes in bp (default: 500 1000)")
    args = parser.parse_args()
    main(run_tests=args.test, window_sizes=args.window_sizes)
