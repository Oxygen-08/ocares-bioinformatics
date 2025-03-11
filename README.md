# Bacterial Genome Analysis Project

## Project Overview
This project provides a set of tools and scripts for analyzing bacterial genomes, focusing on E. coli strains and related species. The analysis includes genome comparison, blast analysis, and identification of divergent regions.

## Project Structure

The project is organized as follows:

```
project_root/
├── data/
│   ├── raw_genomes/          # Original, unmodified genome files
│   ├── processed_genomes/    # Modified or combined genome files
│   └── reference_genomes/    # Reference genome files used for comparisons
├── scripts/
│   ├── analysis/             # Scripts for data analysis
│   ├── visualization/        # Scripts for generating visualizations
│   └── utils/                # Utility scripts for data processing
├── results/
│   ├── alignments/           # Genome alignment files (.coords, .delta)
│   ├── blast_results/        # BLAST comparison outputs
│   ├── divergent_regions/    # Identified divergent regions between genomes
│   └── figures/              # Generated plots and visualizations
├── notebooks/                # Jupyter notebooks for interactive analyses
└── docs/                     # Documentation files
```

## File Naming Conventions

### Genome Files
- Format: `[species]_[serotype/strain]_[additional_info].fasta`
- Example: `ecoli_o157h7_ref.fasta` - Reference genome for E. coli O157:H7

### Analysis Results
- Alignment files: `[reference]_vs_[query].[format]`
- Divergent regions: `[genome]_divergent_regions.[format]`
- BLAST results: `blast_[query]_vs_[reference].[format]`

## Directory Contents

### data/
Contains all genome sequence files used in the project:

- **raw_genomes/**: Original, unmodified genome files for various E. coli strains and related species.
  - Includes strains like O157:H7, O111, O145:H28, O103:H2, and related species like E. fergusonii and E. alberti.

- **processed_genomes/**: Genomes that have been modified or combined for specific analyses.
  - Contains combined_genomes.fasta which includes multiple genomes for comparative analysis.

- **reference_genomes/**: Key reference genomes used as standards for comparisons.
  - Includes E. coli O157:H7, E. coli SE11, and E. coli K12 references.

### scripts/

- **analysis/**:
  - `analyze_blast.py`: Processes BLAST comparison results
  - `blast_statistics.py`: Generates statistical summaries of BLAST comparisons
  - `divergent_analysis.py`: Analyzes divergent regions between genomes
  - `extract_divergent_regions.py`: Extracts sequences from divergent regions

- **visualization/**:
  - `blast_results_heatmap.py`: Generates heatmaps of BLAST comparison results
  - `blast_results_heatmap_adapted.py`: Modified version with additional features

- **utils/**:
  - `count_bases.py`: Utility for counting nucleotide compositions
  - `merge-overlaps.py`: Merges overlapping genomic regions
  - `suffix_tree_findmatch.py`: Finds sequence matches using suffix tree algorithm

### results/

- **alignments/**: Contains genome alignment results
  - Includes .coords files (showing coordinates of aligned regions)
  - .delta files (NUCmer alignment outputs)
  - matches.mum files (MUMmer match coordinates)

- **blast_results/**: BLAST comparison outputs and analyses

- **divergent_regions/**: Contains files identifying regions that differ between genomes
  - BED files marking aligned and non-aligned regions
  - Extracted divergent sequences

- **figures/**: Visualization outputs including heatmaps, phylogenetic trees, and other analysis plots

### notebooks/
Contains Jupyter notebooks (`project.ipynb`) for interactive analysis and result visualization.

### docs/
Contains project documentation, including exported notebook HTML files for easy viewing.

## Getting Started

1. Clone this repository
2. Navigate to the scripts directory to run specific analyses:
   - For BLAST analysis: `python scripts/analysis/analyze_blast.py`
   - For visualization: `python scripts/visualization/blast_results_heatmap.py`
3. Explore the notebooks directory for interactive examples

## Dependencies

- Python 3.6+
- Biopython
- Pandas
- Matplotlib
- Seaborn
- NumPy
- Jupyter (for notebooks)

## Analysis Workflow

1. Genome comparison using NUCmer/MUMmer
2. BLAST analysis to identify sequence similarities
3. Extraction and analysis of divergent regions
4. Visualization of results through heatmaps and other plots

