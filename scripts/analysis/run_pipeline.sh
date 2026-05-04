#!/usr/bin/env bash
# Master pipeline runner — executes all steps sequentially within fp_pipeline env.
# Usage:
#   conda activate fp_pipeline
#   bash scripts/analysis/run_pipeline.sh [--skip-download] [--skip-simulate]
#
# Steps:
#   01  Download E. coli reference genomes
#   02  NUCmer tiered comparative analysis
#   03  Extract discriminatory marker sequences
#   05  Simulate metagenomic community (InSilicoSeq)
#   04  BLAST markers against simulated metagenome
#   06  XGBoost classifier + SHAP explainability
#
# Note: Anvi'o pangenome (step 07) must run in the anvio8 environment — see
#       scripts/analysis/07_pangenome_anvio.sh

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
SCRIPTS="${REPO_ROOT}/scripts/analysis"

# Use explicit env Python to avoid conda run PATH resolution issues
PYTHON="/opt/anaconda3/envs/fp_pipeline/bin/python"
if [[ ! -x "${PYTHON}" ]]; then
    echo "ERROR: fp_pipeline conda environment not found at ${PYTHON}"
    echo "       Run: conda env create -f environment.yml"
    exit 1
fi

SKIP_DOWNLOAD=false
SKIP_SIMULATE=false
for arg in "$@"; do
    [[ "$arg" == "--skip-download"  ]] && SKIP_DOWNLOAD=true
    [[ "$arg" == "--skip-simulate"  ]] && SKIP_SIMULATE=true
done

log() { echo "[$(date +%H:%M:%S)] $*"; }

log "=== Reducing False Positives — Pipeline Run ==="
log "Repository: ${REPO_ROOT}"

if [[ "$SKIP_DOWNLOAD" == false ]]; then
    log "--- Step 01: Download genomes ---"
    "${PYTHON}" "${SCRIPTS}/01_download_genomes.py"
fi

log "--- Step 02: NUCmer tiered analysis ---"
"${PYTHON}" "${SCRIPTS}/02_nucmer_tiered.py"

log "--- Step 03: Extract discriminatory markers ---"
"${PYTHON}" "${SCRIPTS}/03_extract_markers.py"

if [[ "$SKIP_SIMULATE" == false ]]; then
    log "--- Step 05: Simulate metagenome ---"
    "${PYTHON}" "${SCRIPTS}/05_simulate_metagenome.py"
fi

log "--- Step 04: BLAST screening ---"
"${PYTHON}" "${SCRIPTS}/04_blast_screen.py"

log "--- Step 06: ML classifier ---"
"${PYTHON}" "${SCRIPTS}/06_ml_classifier.py"

log "=== Pipeline complete. Results in ${REPO_ROOT}/data/results/ ==="
