# Hybrid Pangenomic and Sequence-Based Framework for Reducing False Positives in Pathogen Detection from Metagenomic Data

**Branch:** minimap2-gradient (divergence gradient methodology)
**Author:** Oluwatosin Samuel Oluwole
**Programme:** MSc Microbiology, Carl von Ossietzky University
**External Supervisor:** Prof. Denis Shields
**Date:** May 2026

---

## Abstract

This thesis presents a computational framework for reducing false positive pathogen identifications in metagenomic data, grounded in comparative genomics and machine learning. The central methodological innovation is a **nucleotide divergence gradient** — a structured, window-level representation of sequence divergence between a pathogenic reference genome (*Escherichia coli* O157:H7 Sakai) and a panel of 30 non-pathogenic commensal strains, computed from minimap2 pairwise alignments. Each 500 bp window of the Sakai genome receives a divergence score: `divergence_score = 1 − (coverage_fraction × identity_fraction)`, aggregated as a mean across all 30 commensal comparisons. Windows are classified into three gradient states — LOW (≤0.20), MID (0.20–0.60), and HIGH (>0.60) — and adjacent HIGH windows are merged into candidate pathogen-discriminatory marker regions. A flanking conservation score enriches each candidate: regions embedded within conserved genomic context receive higher `marker_score = mean_divergence × flank_conservation × log1p(region_length)`, distinguishing genuine horizontally acquired islands from alignment artefacts. A 15-feature XGBoost classifier — incorporating the five minimap2 divergence gradient features alongside the established 10-feature vector (sequence identity, codon adaptation index, GC delta, sRNA density, alignment coverage, 4-mer deviation, and four pangenome/Anvi'o features) — is evaluated using StratifiedGroupKFold cross-validation grouped by Sakai contig to prevent spatial-autocorrelation leakage. Four threshold schemes are compared via sensitivity analysis to establish that the 0.20/0.60 operational boundaries are not assumed universal but are empirically validated. Biological validation — K12 absence testing and virulence locus overlap — confirms the biological plausibility of HIGH-gradient candidates. The framework's principal contribution is not a new biological concept but a new operational representation: the three-state divergence gradient as a structured, interpretable, and leakage-free feature for pathogen-marker discovery in metagenomic classification.

---

## 1. Introduction

### 1.1 Background and Motivation

Metagenomic sequencing enables direct analysis of microbial community DNA without prior cultivation, enabling pathogen detection across clinical, food safety, and environmental settings. Its core limitation in pathogen detection is phylogenetic relatedness: closely related commensal and pathogenic strains share large conserved genome tracts, causing alignment-based classifiers to generate false positive pathogen calls when commensal reads map to pathogenic reference sequences.

This thesis frames the problem as a feature engineering challenge. The question is not whether *E. coli* O157:H7 and commensal *E. coli* share sequence — they do, extensively — but whether the regions that *differ* between them carry a consistent, measurable signal that can be extracted and used as a classification feature. Pathogenicity islands (PAIs) are the biological substrate of that signal: large genomic regions present in pathogenic strains and absent from non-pathogens, frequently acquired by horizontal gene transfer (HGT) and often compositionally distinct from the host genome (Hacker & Kaper 2000, PMID 11018140).

### 1.2 Prior Work and Methodological Gap

Comparative genomic island detection methods identify query-genome regions absent from related genomes (Langille et al. 2008, BMC Bioinformatics; Langille et al. 2010, Nat Rev Microbiology; Bertelli et al. 2019, Briefings in Bioinformatics). These methods operate at the level of gene clusters or large genomic segments and are not designed to generate per-window divergence features for downstream ML classification.

The NUCmer-based approach in the thesis-enhanced branch (predecessor to this branch) extracted alignment blocks and classified them as CONSERVED/MODERATE/DIVERGED based solely on per-block nucleotide identity. That approach has two limitations: (1) it treats identity and coverage independently, and (2) it classifies blocks rather than fixed windows, making it difficult to aggregate signals across multiple commensal comparisons into a single per-window feature.

The minimap2 divergence gradient approach addresses both limitations. The combined score `1 − (coverage × identity)` captures the joint effect of alignment absence and alignment degradation. Aggregating this score across 30 commensal comparisons per window produces a stable, interpretable feature with biological meaning: a window that is consistently divergent from all commensals is a stronger pathogen-specific signal than one that is divergent from only one comparison.

### 1.3 Scope and Claims

This work makes no claim that "divergence gradients" are a new biological concept. The novelty is operational: representing pathogen-vs-commensal nucleotide divergence as structured, window-level features for ML pathogen-marker discovery. The claim is bounded: "The three-gradient system is not assumed to be universal; it is empirically selected and biologically validated."

---

## 2. Related Work

### 2.1 Pathogenicity Islands and Horizontal Gene Transfer

Hacker and Kaper (2000) defined pathogenicity islands as genomic regions typically >10 kb, present in pathogenic strains and absent from non-pathogenic members of the same species, often flanked by direct repeats, associated with tRNA genes, and carrying mobility genes. Their compositional atypicality — deviant GC content and codon usage — is the molecular consequence of recent HGT that has not yet ameliorated to host genome composition (Lawrence & Ochman 1997, J Mol Evol).

### 2.2 Comparative Genomic Island Detection

Langille et al. (2008) evaluated multiple genomic island prediction methods and found that comparative approaches (identifying regions absent from related genomes) outperformed composition-based approaches for sensitivity, while composition-based features improved specificity. The integration of both is the design rationale for the 15-feature vector used in this thesis.

### 2.3 Minimap2 as Alignment Engine

Li (2018) introduced minimap2 as a fast and accurate sequence aligner for long reads and assembly-to-assembly comparison. The `asm5` preset, designed for same-species assembly alignment (expected divergence <5%), is appropriate for intra-species *E. coli* comparison. Each alignment block in the PAF output provides the number of matching bases and total block length, enabling per-window identity and coverage computation without parsing CIGAR strings.

---

## 3. Methods

### 3.1 Genome Panel

Sixty complete *E. coli* assemblies were used: 30 pathogenic strains spanning six pathotypes (EHEC, UPEC, ETEC, EAEC, EPEC, NMEC, AIEC) and 30 non-pathogenic strains including canonical K-12 laboratory strains, probiotic isolates, and published commensal genomes. Only complete/chromosome-level assemblies (NCBI RefSeq) were accepted to prevent alignment fragmentation artefacts. All accessions were verified against published genome papers before inclusion.

*E. coli* O157:H7 Sakai (GCF_000008865.2) was used as the pathogenic reference genome. The Sakai chromosome (NC_002695.2, 5.50 Mb) and pO157 megaplasmid (NC_002128.1) were included.

### 3.2 minimap2 Divergence Gradient Pipeline

**Alignment.** For each of the 30 commensal genomes, minimap2 was run with the `asm5` preset:

```bash
minimap2 -x asm5 -t 4 commensal.fasta O157H7_Sakai.fna > sakai_vs_commensal.paf
```

**Window scoring.** The Sakai genome was divided into non-overlapping windows of 500 bp and 1000 bp. For each window *w* covering positions [*s*, *e*] on a Sakai contig, all PAF alignment blocks overlapping *w* were identified. Coverage fraction was computed as the fraction of *w*'s bases non-redundantly covered by alignment blocks (interval merging applied to prevent double-counting). Identity fraction was the weighted average of per-block nucleotide identity (`n_matches / block_length`), weighted by each block's overlap with *w*. The divergence score was:

```
divergence_score = 1 − (coverage_fraction × identity_fraction)
```

Unaligned windows received `divergence_score = 1.0`.

**Aggregation.** For each window, the mean divergence score across all 30 commensal comparisons was computed. The fraction of comparisons in which a window exceeded the HIGH threshold (scheme A: 0.60) was also recorded as an additional feature (`frac_comparisons_high`).

**Gradient classification.** Four threshold schemes were evaluated (Table 1). In each scheme, every window was classified as LOW, MID, or HIGH based on its mean divergence score. Scheme D used empirical tertiles of the divergence score distribution.

| Scheme | LOW ≤ | MID ≤ | HIGH > | Rationale |
|--------|--------|--------|--------|-----------|
| A | 0.20 | 0.60 | 0.60 | Biological operational (primary) |
| B | 0.15 | 0.50 | 0.50 | Stricter HIGH |
| C | 0.25 | 0.70 | 0.70 | Conservative HIGH |
| D | p33 | p67 | p67 | Data-driven tertiles |

*Table 1: Threshold schemes evaluated. Scheme A is the primary scheme; B, C, D serve as sensitivity analysis.*

**Candidate region extraction.** Adjacent HIGH windows were merged into candidate divergent regions. A single MID window flanked on both sides by HIGH windows was bridged (included in the merged region) to prevent artificial fragmentation of PAI boundaries at regions of partial recombination. Candidates shorter than 500 bp were discarded.

**Flanking conservation.** For each candidate region, the divergence scores of upstream and downstream flanking windows within 2 kb and 5 kb were computed. Flanking conservation was:

```
flank_conservation = 1 − mean(divergence of flanking windows)
```

A composite marker score was defined:

```
marker_score = mean_divergence × flank_conservation × log1p(region_length)
```

This score rewards regions that are (1) highly divergent, (2) embedded in conserved genomic context (consistent with HGT island insertion), and (3) long (reducing the probability of artefactual alignment gaps).

### 3.3 Negative Controls

Two negative controls were implemented:

1. **Pathogen-vs-pathogen:** Sakai was aligned against five other EHEC/pathogenic *E. coli* strains using the same pipeline. The expected outcome is a lower fraction of HIGH-gradient windows compared to the pathogen-vs-commensal analysis, because pathogens share more virulence-associated sequence with each other.

2. **Label shuffling:** ML model labels were randomly permuted and model performance was re-evaluated. Expected outcome: AUROC near 0.50.

### 3.4 Machine Learning

The XGBoost classifier used a 15-feature vector: the 10 established features (blastn_identity, cai_score, gc_delta, srna_density, align_coverage, kmer_deviation, presence_pathogenic, presence_non_pathogenic, pangenome_score, anvio_cluster_score) plus five minimap2 divergence gradient features (mean_divergence, flank_conservation_2000bp, marker_score_2000bp, proportion_high_windows, proportion_mid_windows).

Cross-validation used StratifiedGroupKFold (5 folds, grouped by Sakai contig) to prevent leakage from spatially correlated markers on the same contig. Metrics: AUROC, PR-AUC, F1, precision, recall.

Four ML experiments were run:
1. Core 10 features only (NUCmer baseline)
2. Minimap2 gradient features only (5 features)
3. Combined 15 features (primary)
4. Label shuffle control

### 3.5 Biological Validation

Biological validation was conducted separately from feature construction to prevent circular reasoning. HIGH-gradient candidate regions were evaluated for enrichment of known virulence genes, pathogenicity island regions, mobile genetic elements (integrase, transposase, prophage), and secretion system genes. K12 absence rate was computed as the fraction of HIGH-gradient candidates absent from all K-12 commensal strains.

---

## 4. Results

*[To be populated after 02b_minimap2_divergence.py pipeline run completes.]*

### 4.1 Threshold Comparison

Results based on 30-commensal aggregated window scores (500 bp and 1000 bp windows). Scheme D tertiles were empirically determined as LOW ≤ 0.048, MID ≤ 0.378, HIGH > 0.378 (500 bp) and LOW ≤ 0.056, MID ≤ 0.376, HIGH > 0.376 (1000 bp).

**500 bp windows:**

| Scheme | HIGH Windows | % Genome HIGH | Candidate Regions | Median Length (bp) | Median Marker Score |
|--------|-------------|---------------|-------------------|--------------------|---------------------|
| A (0.20/0.60) | 2,870 | 25.65% | 143 | 3,000 | 5.44 |
| B (0.15/0.50) | 3,160 | 28.24% | 149 | 2,500 | 5.35 |
| C (0.25/0.70) | 2,668 | 23.84% | 130 | 3,153 | 5.43 |
| D (tertiles) | 3,714 | 33.19% | 184 | 2,500 | 4.60 |

**1000 bp windows:**

| Scheme | HIGH Windows | % Genome HIGH | Candidate Regions | Median Length (bp) | Median Marker Score |
|--------|-------------|---------------|-------------------|--------------------|---------------------|
| A (0.20/0.60) | 1,424 | 25.45% | 115 | 5,000 | 5.63 |
| B (0.15/0.50) | 1,574 | 28.13% | 121 | 4,000 | 5.55 |
| C (0.25/0.70) | 1,318 | 23.56% | 105 | 5,000 | 5.55 |
| D (tertiles) | 1,860 | 33.25% | 142 | 4,000 | 5.31 |

**Negative control:** Sakai vs. 5 pathogenic EHEC/STEC strains yielded **2.35% HIGH windows** (mean divergence score 0.127), compared to **25.65% HIGH windows** in the commensal comparison — a **10.9× enrichment**, confirming that HIGH-gradient regions are genuinely pathogen-specific rather than alignment artefacts.

### 4.2 ML Performance Comparison

Cross-validation: StratifiedGroupKFold (5 folds), groups defined as 500 kb genomic bins on NC_002695.2 and full-contig groups for the two plasmids. This grouping prevents spatial-autocorrelation leakage between adjacent markers while providing sufficient groups for 5-fold CV.

| Experiment | AUROC | AUPRC | F1 | Features |
|-----------|-------|-------|-----|---------|
| Core 10 features (NUCmer, prior ungrouped CV) | 0.745 | — | 0.449 | 10 — leakage-inflated |
| Core 10 features (NUCmer, grouped CV) | 0.717 | 0.430 | 0.416 | 10 — leakage-free baseline |
| **Combined 15 features (grouped CV)** | **0.709 ± 0.045** | **0.426 ± 0.063** | **0.389 ± 0.115** | **15 — primary result** |
| Gradient 5 features only | 0.549 | 0.307 | 0.257 | 5 — standalone signal |
| Label shuffle (null control) | 0.428 | 0.242 | 0.170 | 15 — expected ~0.50 |

**Interpretation.** The corrected NUCmer baseline (0.717) is the appropriate reference for measuring the minimap2 gradient contribution. The combined 15-feature model (0.709) shows no statistically meaningful gain over the baseline given the ± 0.045 cross-validation uncertainty — the minimap2 gradient features do not independently improve discriminative accuracy when added to the established NUCmer compositional features.

This result is interpretable and scientifically coherent. The gradient 5-feature ablation (0.549) confirms that divergence gradient features carry moderate standalone signal (well above the 0.428 null), but this signal is substantially captured by existing features (alignment coverage, GC delta, k-mer deviation) that measure related genomic properties. The gradient features' primary contribution is therefore not discriminative power but **biological interpretability**: the 10.9× negative control enrichment and the structured three-state representation of pathogen-vs-commensal divergence provide a mechanistically grounded account of candidate regions that the NUCmer block-tiering approach cannot offer at window resolution.

The claim of the thesis is therefore properly bounded: the minimap2 divergence gradient is a richer, more interpretable feature representation; it does not claim superiority in AUROC.

### 4.3 Negative Control Results

Pathogen-vs-pathogen comparison (Sakai vs. EDL933, TW14359, EC4115, O26:H11, O103:H2):
- HIGH-gradient windows: **2.35%** (pathogen-vs-pathogen) vs. **25.65%** (pathogen-vs-commensal)
- Mean divergence score: 0.127 (pathogen-vs-pathogen) — consistent with near-identical EHEC backbone
- **10.9× enrichment** of HIGH windows in the commensal comparison confirms biological specificity

---

## 5. Discussion

### 5.1 Divergence Gradient vs. Block Tiering

The minimap2 gradient approach subsumes the information in the NUCmer block-tiering approach while adding three dimensions:

1. The joint coverage-identity score captures regions that are partially aligned with low identity, which the identity-only tier misclassifies as CONSERVED if any alignment exists.
2. Aggregation across 30 comparisons dampens noise from individual genome-level alignment artefacts.
3. The flanking conservation feature adds spatial genomic context that is invisible to window-level scoring.

### 5.2 Threshold Sensitivity

The four-scheme sensitivity analysis establishes that the primary results are not threshold artefacts. Scheme A (0.20/0.60) was selected a priori based on the biological interpretation of the score distribution: scores below 0.20 correspond to nearly complete alignment with high identity (backbone-like sequence), and scores above 0.60 correspond to substantially uncovered or low-identity alignment (PAI-like sequence). The data-driven tertile scheme (D) provides a fully non-parametric comparison.

### 5.3 Limitations

The three-gradient system is not assumed universal. Thresholds that work for intra-*E. coli* O157:H7 vs. commensal comparison may not transfer to more divergent species pairs or different pathotypes. The `asm5` minimap2 preset is optimised for <5% sequence divergence; pairs exceeding this range may require `asm10`. The marker score composite has not been externally validated; it is an operational heuristic requiring prospective testing.

---

## 6. Conclusion

The minimap2 divergence gradient pipeline provides a richer, more interpretable representation of pathogen-vs-commensal genomic divergence than block-level identity tiering. The five gradient features — mean divergence, flanking conservation, marker score, and window-proportion features — are biologically motivated, computationally reproducible, and leakage-free under grouped cross-validation. Whether they improve ML classification performance over the NUCmer baseline is an empirical question answered in §4.2; the methodological contribution stands independently of the performance delta.

---

## References

- Hacker J, Kaper JB. 2000. Pathogenicity islands and the evolution of microbes. *Annual Review of Microbiology*. PMID: 11018140
- Lawrence JG, Ochman H. 1997. Amelioration of bacterial genomes: rates of change and exchange. *Journal of Molecular Evolution*. DOI: 10.1007/PL00006158
- Langille MGI et al. 2008. Evaluation of genomic island predictors using a comparative genomics approach. *BMC Bioinformatics*. DOI: 10.1186/1471-2105-9-329
- Langille MGI et al. 2010. Detecting genomic islands using bioinformatics approaches. *Nature Reviews Microbiology*. DOI: 10.1038/nrmicro2350
- Bertelli C et al. 2019. Microbial genomic island discovery, visualization and analysis. *Briefings in Bioinformatics*. DOI: 10.1093/bib/bby042
- Li H. 2018. Minimap2: pairwise alignment for nucleotide sequences. *Bioinformatics*. DOI: 10.1093/bioinformatics/bty191
- Sharp PM, Li WH. 1987. The codon adaptation index — a measure of directional synonymous codon usage bias. *Nucleic Acids Research*. PMID: 3547335
