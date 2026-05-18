# Hybrid Pangenomic Framework for Reducing False Positives in Metagenomic Pathogen Detection

**Author:** Oluwatosin Samuel Oluwole  
**Programme:** MSc Microbiology, Carl von Ossietzky University  
**Supervisor:** Prof. Denis Shields  

Computational pipeline for the MSc thesis. Identifies *E. coli* O157:H7-specific pathogenicity markers through NUCmer tiered comparative genomics across a balanced 30/30 pathogenic/non-pathogenic panel, validates them population-wide via Anvi'o pangenomics, and classifies them using a 15-feature XGBoost model (AUROC 0.724 ± 0.028 vs BLAST baseline 0.552).

---

## Repository Structure

```
bioinformatics/
├── scripts/analysis/           # Full pipeline — numbered in execution order
│   ├── 01_download_genomes.py  # Download assemblies from NCBI via datasets CLI
│   ├── 02_nucmer_tiered.py     # NUCmer pairwise alignment + tier classification
│   ├── 02b_minimap2_divergence.py  # minimap2 divergence gradient (Layer 6 features)
│   ├── 03_extract_markers.py   # Extract 415 candidate markers from alignment blocks
│   ├── 04_simulate_metagenome.py   # InSilicoSeq simulated metagenome (HiSeq model)
│   ├── 05_blast_screen.py      # BLASTn screen + AUROC baseline
│   ├── 06_pangenome_anvio.sh   # Anvi'o pangenome — separate env, run BEFORE step 07
│   ├── 07_ml_classifier.py     # 15-feature XGBoost + grouped 5-fold CV
│   ├── 08_bio_sanity.py        # K-12 absence test + virulence locus overlap
│   ├── 08b_divergence_landscape.py # Chromosome-wide divergence landscape figures
│   ├── 09_thesis_slides.py     # Auto-generate thesis presentation (15 slides)
│   ├── 10_kraken2_comparison.py    # Kraken2 custom-DB benchmark
│   ├── 11_triangulate_experiments.py  # Regularisation, K-12 relabelling, null control
│   ├── fig_alignment_landscape.py  # Figure: alignment landscape
│   ├── fig_pangenome.py            # Figure: pangenome composition + enrichment
│   ├── fig_results.py              # Figure: ROC/PR curves, SHAP, confusion matrices
│   ├── fix_read_labels.py          # Utility: fix InSilicoSeq read label bug
│   ├── run_missing_nucmer.py       # Utility: append NUCmer blocks for new commensal strains
│   ├── run_new_pathogens_nucmer.py # Utility: append NUCmer blocks for new pathogenic strains
│   └── run_pipeline.sh             # Orchestrator: steps 01–08 + figures
├── data/
│   ├── genomes/
│   │   └── genome_manifest.tsv # 60-strain manifest (GCF accessions, pathotypes, PMIDs)
│   ├── metagenome/             # Simulated reads (gitignored — regenerate via step 04)
│   └── results/                # Key TSV outputs and figures (tracked)
├── thesis/
│   ├── thesis_v2_synthesis.md       # FINAL SUBMISSION — MSc thesis
│   ├── project1_report.md           # Report 1: Comparative Genomics (Short Project 1)
│   ├── project2_report.md           # Report 2: Pangenomics (Short Project 2)
│   ├── thesis_minimap2_gradient.md  # Exploration: minimap2 as standalone approach (not submitted)
│   ├── thesis_v1_draft.md           # Historical: first draft (superseded)
│   └── pathogen-detection-workflow.drawio  # Pipeline diagram
├── notebooks/
│   └── 00_prototype_exploration.ipynb
├── environment.yml             # fp_pipeline conda environment (main pipeline)
├── environment_anvio.yml       # Anvi'o conda environment (step 06 only)
└── .pre-commit-config.yaml     # nbstripout hook
```

---

## Reproducing the Pipeline

### Prerequisites

Two conda environments are required — the main pipeline and Anvi'o (strict dependency isolation):

```bash
# Main pipeline
conda env create -f environment.yml
conda activate fp_pipeline

# Anvi'o (step 06 only)
conda env create -f environment_anvio.yml
conda activate anvio8
```

### Execution Order

Step 06 (Anvi'o pangenome) must complete **before** Step 07 (ML classifier), because
`07_ml_classifier.py` reads `anvio_cluster_score` from the Anvi'o output.

**Steps 01–05, 07–11 and all figures** (fp_pipeline env):

```bash
conda activate fp_pipeline
bash scripts/analysis/run_pipeline.sh
```

Options:
```
--skip-download    Skip genome download if assemblies already present
--skip-simulate    Skip InSilicoSeq simulation if reads already exist
--skip-plots       Skip figure generation
```

**Step 06 — Anvi'o Pangenome** (separate environment, run first):

```bash
conda activate anvio8
bash scripts/analysis/06_pangenome_anvio.sh
```

This step requires ~4–6 hours on a modern laptop. Outputs are written to `data/results/pangenome/`.

### Downloading Genomes

Genomes are not stored in the repository (gitignored). Re-download using the manifest:

```bash
conda activate fp_pipeline
python scripts/analysis/01_download_genomes.py
```

All 60 GCF accessions are recorded in `data/genomes/genome_manifest.tsv` with pathotypes and supporting PMIDs.

---

## Data Availability

Large files (genome assemblies, simulated reads, NUCmer delta files) are gitignored. Key results tracked in the repository:

| File | Description |
|------|-------------|
| `data/genomes/genome_manifest.tsv` | 60-strain manifest: GCF accessions, pathotypes, PMIDs, FASTA paths |
| `data/results/tiered_blocks.tsv` | 23,923 NUCmer alignment blocks (60 strains vs Sakai) |
| `data/results/ml/feature_matrix.tsv` | 415 markers × 15 features |
| `data/results/ml/cv_results.tsv` | 5-fold grouped CV metrics (per fold) |
| `data/results/ml/shap_values.tsv` | Per-feature SHAP values (415 markers) |
| `data/results/blast_screen/classification_metrics.tsv` | BLAST baseline per tier |
| `data/results/pangenome/enrichment_scores.tsv` | COG14 enrichment (2,620 rows) |
| `data/results/tier_summary.tsv` | 415-marker tier breakdown |
| `data/results/figures/*.png` | All thesis figures |

---

## Key Results

| Metric | Value |
|--------|-------|
| Comparison panel | 60 strains (30 pathogenic / 30 non-pathogenic, balanced) |
| NUCmer alignment blocks | 23,923 |
| Candidate markers extracted | 415 (134 CONSERVED / 172 MODERATE / 109 DIVERGED) |
| BLAST screen AUROC (baseline) | 0.552 |
| Kraken2 custom-DB AUROC (baseline) | 0.513 |
| ML classifier AUROC (15-feature XGBoost, grouped CV) | **0.724 ± 0.028** |
| ML classifier AUPRC | **0.450 ± 0.041** (random baseline = 0.263) |
| DIVERGED markers absent from all K-12 strains | 54 / 109 (49.5%) |
| Markers overlapping named EHEC virulence loci | 14 (LEE, Stx1/2, OI-48, tellurite) |
| Anvi'o pangenome gene clusters | 12,204 (46-genome Anvi'o panel) |
| COG14 functions enriched in pathogens (q < 0.05) | 6 (T3SS, phage, transposases) |
| minimap2 negative control enrichment | 10.9× (HIGH windows: commensal vs pathogen-vs-pathogen) |

---

## Citation

Oluwole, O. S. (2026). *Hybrid Pangenomic and Sequence-Based Framework for Reducing False Positives in Pathogen Detection from Metagenomic Data*. MSc Thesis, Carl von Ossietzky University.
