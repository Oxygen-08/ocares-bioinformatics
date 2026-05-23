# MASTER DIRECTIVE: Dual-Feature Evolutionary Integration Pipeline (RSCU Cosine Distance & Delta-CAI)

## 🧬 Context & Background
We are building a machine learning classifier to detect E. coli O157:H7 metagenomic signals. Our dataset includes a 61-genome pangenome yielding 15,251 total gene clusters (1,962 core, 6,189 accessory). 

To optimize classifier performance (current AUROC 0.711) and provide explainable biological data for our thesis text, we are implementing two distinct features based on evolutionary amelioration theory:
1. **Cosine Distance of RSCU Vectors (ML Feature):** Captures multi-dimensional codon style for high-density ML gradients.
2. **Delta-CAI (Biological Feature):** Captures absolute translation velocity differences for academic discussion.

## 🛠 Target Architecture & Steps
You must work inside the active Conda environment (`fp_pipeline`) using Biopython. Execute these phases sequentially. Do not skip testing. Commit files to git as you complete each phase.

### Phase 1: Core Genome Profiling & Baseline Matrices
1. Write a script: `scripts/analysis/calculate_core_reference.py`.
2. This script must read the nucleotide sequences of the 1,962 CORE gene clusters (extracted from our Anvi'o pangenome).
3. Generate two outputs from this core baseline:
   - A single **CAI Reference Weight Matrix** (based on Sharp & Li, 1987).
   - A master **59-dimensional Core RSCU Reference Vector** (excluding stop codons, Methionine [ATG], and Tryptophan [TGG]).
4. Save these profiles as a single unified JSON file in `data/metrics/core_evolutionary_reference.json`.

### Phase 2: Dual-Feature Engineering Script
1. Write a script: `scripts/analysis/compute_evolutionary_features.py`.
2. Load the JSON baseline reference file from Phase 1.
3. For every gene sequence within the 6,189 accessory clusters (prioritizing the 66 pathogen-enriched COG clusters):
   - **Calculate the CAI score** and subtract the median core CAI to get `Delta_CAI` (Scalar).
   - **Calculate the 59-dimensional RSCU vector** for that gene, and compute the mathematical `Cosine Distance` between it and the core vector.
4. Output these two metrics for every accessory gene cluster.

### Phase 3: Matrix Integration & Dataframe Update
1. Update our master ML feature matrix (`data/results/ml/feature_matrix.tsv`).
2. Add `Cosine_Distance` and `Delta_CAI` as new continuous columns for every gene present in a genome.
3. For missing/absent genes within a genome profile, input a standardized neutral flag value (e.g., maximum distance or `1.0`) and document your choice clearly in the script comments.

### Phase 4: Statistical Validation & Reporting
1. Print out a summary statistic showing:
   - Variance and mean of `Cosine_Distance` in pathogen-enriched clusters vs. commensal-enriched/remaining accessory clusters.
   - Variance and mean of `Delta_CAI` across the same splits.
