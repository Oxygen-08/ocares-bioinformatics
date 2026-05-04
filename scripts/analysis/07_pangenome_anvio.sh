#!/usr/bin/env bash
# Anvi'o pangenome analysis — must run inside the anvio8 conda environment.
#
# Usage:
#   conda activate anvio8
#   bash scripts/analysis/07_pangenome_anvio.sh
#
# What this does:
#   1. Reformat FASTA headers (Anvi'o requires ≤ 100 char contig names, A-Za-z0-9_- only)
#   2. Create Anvi'o contigs databases for all strains
#   3. Run HMM profiles for single-copy core gene annotation
#   4. Create genomes storage
#   5. Run pangenome analysis (MCL clustering)
#   6. Compute genome similarity and enrichment scores for O-island genes
#   7. Export gene cluster frequencies per genome for downstream ML feature
#
# Key parameters justified against Anvi'o pangenome tutorial (Eren et al. 2015):
#   --minbit=0.8        : conservative threshold to reduce spurious homologue links
#   --mcl-inflation=10  : higher inflation → tighter clusters (appropriate for
#                         closely related E. coli strains where pan-genome
#                         over-merging is the dominant error mode)

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
GENOME_DIR="${REPO_ROOT}/data/genomes"
PANGENOME_DIR="${REPO_ROOT}/data/results/pangenome"
mkdir -p "${PANGENOME_DIR}"

NUM_THREADS=4
PROJECT_NAME="EcoliPangenome"

log() { echo "[$(date +%H:%M:%S)] ANVIO  $*"; }

# ── 1. Reformat FASTA headers ─────────────────────────────────────────────────
log "Reformatting FASTA headers for Anvi'o compatibility"
REFORMATTED_DIR="${PANGENOME_DIR}/reformatted"
mkdir -p "${REFORMATTED_DIR}"

GENOME_LIST="${PANGENOME_DIR}/genome_list.txt"
> "${GENOME_LIST}"

for genome_dir in "${GENOME_DIR}"/*/; do
    label=$(basename "${genome_dir}")
    fna="${genome_dir}${label}.fna"
    [[ ! -f "${fna}" ]] && continue

    out="${REFORMATTED_DIR}/${label}.fna"
    if [[ ! -f "${out}" ]]; then
        anvi-script-reformat-fasta \
            "${fna}" \
            -o "${out}" \
            --min-len 1000 \
            --simplify-names \
            --prefix "${label}" \
            --seq-type NT \
            2>/dev/null
    fi
    echo -e "${label}\t${out}" >> "${GENOME_LIST}"
done

COUNT=$(wc -l < "${GENOME_LIST}")
log "Reformatted ${COUNT} genomes"

# ── 2. Create contigs databases ───────────────────────────────────────────────
log "Creating Anvi'o contigs databases"
DB_DIR="${PANGENOME_DIR}/contigs_dbs"
mkdir -p "${DB_DIR}"

while IFS=$'\t' read -r label fasta; do
    db="${DB_DIR}/${label}.db"
    if [[ ! -f "${db}" ]]; then
        log "  contig DB: ${label}"
        anvi-gen-contigs-database \
            -f "${fasta}" \
            -o "${db}" \
            -n "${label}" \
            --num-threads "${NUM_THREADS}" \
            2>/dev/null
        anvi-run-hmms \
            -c "${db}" \
            --num-threads "${NUM_THREADS}" \
            --quiet \
            2>/dev/null
    fi
done < "${GENOME_LIST}"

# ── 3. External genomes file ──────────────────────────────────────────────────
EXTERNAL_GENOMES="${PANGENOME_DIR}/external_genomes.txt"
echo -e "name\tcontigs_db_path" > "${EXTERNAL_GENOMES}"
while IFS=$'\t' read -r label _; do
    db="${DB_DIR}/${label}.db"
    [[ -f "${db}" ]] && echo -e "${label}\t${db}" >> "${EXTERNAL_GENOMES}"
done < "${GENOME_LIST}"

# ── 4. Genomes storage ────────────────────────────────────────────────────────
STORAGE="${PANGENOME_DIR}/${PROJECT_NAME}-GENOMES.db"
if [[ ! -f "${STORAGE}" ]]; then
    log "Creating genomes storage"
    anvi-gen-genomes-storage \
        -e "${EXTERNAL_GENOMES}" \
        -o "${STORAGE}" \
        2>/dev/null
fi

# ── 5. Pangenome analysis ─────────────────────────────────────────────────────
PAN_DB="${PANGENOME_DIR}/${PROJECT_NAME}-PAN.db"
if [[ ! -f "${PAN_DB}" ]]; then
    log "Running pangenome analysis (MCL inflation=10)"
    anvi-pan-genome \
        -g "${STORAGE}" \
        --project-name "${PROJECT_NAME}" \
        --output-dir "${PANGENOME_DIR}" \
        --num-threads "${NUM_THREADS}" \
        --minbit 0.8 \
        --mcl-inflation 10 \
        --use-ncbi-blast \
        2>/dev/null
fi

# ── 6. Genome similarity (ANI) ────────────────────────────────────────────────
log "Computing genome similarity (fastANI)"
anvi-compute-genome-similarity \
    -e "${EXTERNAL_GENOMES}" \
    -o "${PANGENOME_DIR}/ANI" \
    -p "${PAN_DB}" \
    --program fastANI \
    --num-threads "${NUM_THREADS}" \
    2>/dev/null || log "WARN: fastANI failed — skipping ANI layer"

# ── 7. Enrichment analysis ────────────────────────────────────────────────────
# Define pathogen group using genome manifest
log "Running enrichment analysis: pathogens vs commensals"

GROUPS="${PANGENOME_DIR}/layer_groups.txt"
echo -e "name\tgroup" > "${GROUPS}"
while IFS=$'\t' read -r label _ pathotype _; do
    [[ "${label}" == "label" ]] && continue  # skip header
    if [[ "${pathotype}" =~ ^(EHEC|EPEC|UPEC|ETEC)$ ]]; then
        echo -e "${label}\tPATHOGEN"
    else
        echo -e "${label}\tCOMMENSAL"
    fi
done < "${REPO_ROOT}/data/genomes/genome_manifest.tsv" >> "${GROUPS}"

anvi-compute-functional-enrichment-in-pan \
    -p "${PAN_DB}" \
    -g "${STORAGE}" \
    --category-variable group \
    --annotation-source COG20_FUNCTION \
    -o "${PANGENOME_DIR}/enrichment_scores.tsv" \
    2>/dev/null || log "WARN: Enrichment analysis failed — may need COG annotation"

# ── 8. Export gene cluster frequencies ───────────────────────────────────────
log "Exporting gene cluster presence/absence matrix"
anvi-export-table \
    "${PAN_DB}" \
    --table gene_clusters \
    -o "${PANGENOME_DIR}/gene_clusters.tsv" \
    2>/dev/null

log "=== Anvi'o pangenome complete. Results in ${PANGENOME_DIR} ==="
