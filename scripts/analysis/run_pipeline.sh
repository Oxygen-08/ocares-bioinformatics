#!/usr/bin/env bash
# ─────────────────────────────────────────────────────────────────────────────
# Master pipeline runner — executes all Python steps sequentially.
#
# Usage:
#   conda activate fp_pipeline
#   bash scripts/analysis/run_pipeline.sh [options]
#
# Options:
#   --skip-download    Skip Step 01 (genome download) if genomes already present
#   --skip-simulate    Skip Step 04 (InSilicoSeq simulation) if reads exist
#   --skip-plots       Skip figure generation (Steps 08, plot_*)
#
# Full pipeline steps:
#   01  Download E. coli reference genomes from NCBI RefSeq
#   02  NUCmer tiered alignment (CONSERVED / MODERATE / DIVERGED)
#   03  Extract discriminatory marker sequences (415 markers)
#   04  Simulate metagenomic community (InSilicoSeq 2.0)
#   05  BLAST screen markers against simulated metagenome
#   07  XGBoost ML classifier + SHAP explainability
#   08  Biological sanity checks (K-12 absence, virulence locus overlap)
#   P1  Alignment landscape figure
#   P2  ML result figures (ROC, confusion matrix, feature correlation)
#   P3  Pangenome figures (composition, enrichment, heatmap)
#
# NOTE: Step 06 (Anvi'o pangenome) runs in a separate conda environment
#       and MUST be completed before Step 07 (ML classifier), because
#       07_ml_classifier.py reads anvio_cluster_score from Anvi'o output.
#       See scripts/analysis/06_pangenome_anvio.sh and README.md.
# ─────────────────────────────────────────────────────────────────────────────
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
SCRIPTS="${REPO_ROOT}/scripts/analysis"

PYTHON="/opt/anaconda3/envs/fp_pipeline/bin/python"
if [[ ! -x "${PYTHON}" ]]; then
    echo "ERROR: fp_pipeline conda environment not found."
    echo "       Run: conda env create -f environment.yml"
    exit 1
fi

SKIP_DOWNLOAD=false
SKIP_SIMULATE=false
SKIP_PLOTS=false
for arg in "$@"; do
    [[ "$arg" == "--skip-download" ]] && SKIP_DOWNLOAD=true
    [[ "$arg" == "--skip-simulate" ]] && SKIP_SIMULATE=true
    [[ "$arg" == "--skip-plots"    ]] && SKIP_PLOTS=true
done

log() { echo "[$(date +%H:%M:%S)] $*"; }

log "=== E. coli O157:H7 False-Positive Reduction Pipeline ==="
log "Repository: ${REPO_ROOT}"
log ""

$SKIP_DOWNLOAD || { log "Step 01 — Download genomes"; "${PYTHON}" "${SCRIPTS}/01_download_genomes.py"; }

log "Step 02 — NUCmer tiered alignment"
"${PYTHON}" "${SCRIPTS}/02_nucmer_tiered.py"

log "Step 03 — Extract discriminatory markers"
"${PYTHON}" "${SCRIPTS}/03_extract_markers.py"

$SKIP_SIMULATE || { log "Step 04 — Simulate metagenome (InSilicoSeq)"; "${PYTHON}" "${SCRIPTS}/04_simulate_metagenome.py"; }

log "Step 05 — BLAST screen"
"${PYTHON}" "${SCRIPTS}/05_blast_screen.py"

log "Step 07 — ML classifier (XGBoost + SHAP)"
"${PYTHON}" "${SCRIPTS}/07_ml_classifier.py"

log "Step 08 — Biological sanity checks"
"${PYTHON}" "${SCRIPTS}/08_bio_sanity.py"

if [[ "$SKIP_PLOTS" == false ]]; then
    log "Plots — Alignment landscape"
    "${PYTHON}" "${SCRIPTS}/plot_alignment_landscape.py"
    log "Plots — ML results (ROC, confusion matrix, feature correlation)"
    "${PYTHON}" "${SCRIPTS}/plot_results.py"
    log "Plots — Pangenome (composition, enrichment, heatmap)"
    "${PYTHON}" "${SCRIPTS}/plot_pangenome.py"
fi

log ""
log "=== Pipeline complete. Results: ${REPO_ROOT}/data/results/ ==="
log "    Step 06 (Anvi'o pangenome): run separately in anvio8 env BEFORE Step 07 — see README.md"
