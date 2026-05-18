# Hybrid Pangenomic Framework for Reducing False Positives in Metagenomic Pathogen Detection

**Author:** Oluwatosin Samuel Oluwole  
**Programme:** MSc Microbiology, Carl von Ossietzky University  
**Supervisor:** Prof. Denis Shields  

Computational pipeline for the master's thesis. Identifies *E. coli* O157:H7-specific pathogenicity markers through NUCmer tiered comparative genomics, validates them population-wide via Anvi'o pangenomics, and classifies them using a 10-feature XGBoost model (AUROC 0.673 vs BLAST baseline 0.552).

---

## Repository Structure

```
bioinformatics/
├── scripts/analysis/          # Full pipeline — run in order
│   ├── 01_download_genomes.py
│   ├── 02_nucmer_tiered.py
│   ├── 03_extract_markers.py
│   ├── 04_simulate_metagenome.py
│   ├── 05_blast_screen.py
│   ├── 06_pangenome_anvio.sh   # Anvi'o only — separate env, run BEFORE step 07
│   ├── 07_ml_classifier.py     # Requires Anvi'o output from step 06
│   ├── 08_bio_sanity.py
│   ├── plot_alignment_landscape.py
│   ├── plot_results.py
│   ├── plot_pangenome.py
│   └── run_pipeline.sh        # Orchestrates steps 01–08 + plots
├── data/
│   ├── genomes/               # Downloaded assemblies (gitignored)
│   ├── metagenome/            # Simulated reads (gitignored)
│   └── results/               # TSV outputs and figures (key files tracked)
├── thesis/
│   ├── thesis_v2_synthesis.md        # Full MSc thesis
│   ├── report_1_comparative_genomics.md
│   ├── report_2_pangenome.md
│   └── thesis_corrected.md           # Original draft (preserved)
├── notebooks/
│   └── project.ipynb
├── mcp_server/                # NCBI E-utilities MCP server (dev tool)
│   └── ncbi_server.py
├── environment.yml            # fp_pipeline conda environment
└── .pre-commit-config.yaml    # nbstripout hook
```

---

## Reproducing the Pipeline

### Prerequisites

Two conda environments are required — the main pipeline and Anvi'o (which has strict dependency isolation):

```bash
# 1. Main pipeline environment
conda env create -f environment.yml
conda activate fp_pipeline

# 2. Anvi'o environment (Step 06 only) — pinned versions
conda env create -f environment_anvio.yml
```

### Execution order

Step 06 (Anvi'o pangenome) must complete **before** Step 07 (ML classifier), because
`07_ml_classifier.py` reads `anvio_cluster_score` from the Anvi'o pangenome output.

**1. Run steps 01–05, 07–08 and all figures** (fp_pipeline env):

```bash
conda activate fp_pipeline
bash scripts/analysis/run_pipeline.sh
```

Options:
```
--skip-download    Skip genome download if already present
--skip-simulate    Skip InSilicoSeq simulation if reads already exist
--skip-plots       Skip figure generation
```

**2. Run Step 06 — Anvi'o Pangenome** (separate environment, run first):

```bash
conda activate anvio8
bash scripts/analysis/06_pangenome_anvio.sh
```

This step requires ~4–6 hours on a modern laptop. Outputs are written to `data/results/pangenome/`.

### Data Availability

Large files (genomes, simulated reads, alignment outputs) are gitignored. Key results tracked in the repo:

| File | Description |
|------|-------------|
| `data/genomes/genome_manifest.tsv` | 60-strain manifest with GCF accessions, pathotypes, PMIDs |
| `data/results/ml/feature_matrix.tsv` | 415 markers × 10 features |
| `data/results/ml/cv_results.tsv` | 5-fold CV metrics |
| `data/results/ml/shap_values.tsv` | Per-feature SHAP values |
| `data/results/blast_screen/classification_metrics.tsv` | BLAST baseline per tier |
| `data/results/pangenome/enrichment_scores.tsv` | COG14 enrichment (2,620 rows) |
| `data/results/tier_summary.tsv` | 415-marker tier breakdown |
| `data/results/figures/*.png` | All thesis figures |

---

## Key Results

| Metric | Value |
|--------|-------|
| Alignment blocks (45 strains vs Sakai) | 17,519 |
| Candidate markers extracted | 415 (134 CONSERVED / 172 MODERATE / 109 DIVERGED) |
| BLAST screen AUROC (baseline) | 0.552 |
| ML classifier AUROC (10-feature XGBoost) | **0.673 ± 0.063** |
| DIVERGED markers absent from all K-12 strains | 54 / 109 (49.5%) |
| Markers overlapping named EHEC virulence loci | 14 (LEE, Stx1/2, OI-48, tellurite) |
| Anvi'o pangenome gene clusters | 12,204 (46 genomes) |
| PATHOGEN-enriched COG14 functions (q < 0.05) | 6 |

---

## NCBI MCP Server (dev tool)

`mcp_server/ncbi_server.py` wraps NCBI E-utilities for live database queries within Claude Code sessions. Not required to reproduce the pipeline.

```bash
export NCBI_EMAIL=your@email.com
pip install -r mcp_server/requirements.txt
```

---

## Citation

Oluwole, O. S. (2025). *Hybrid Pangenomic and Sequence-Based Framework for Reducing False Positives in Pathogen Detection from Metagenomic Data*. MSc Thesis, Carl von Ossietzky University.
