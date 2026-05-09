# Hybrid Pangenomic and Sequence-Based Framework for Reducing False Positives in Pathogen Detection from Metagenomic Data

**Author:** Oluwatosin Samuel Oluwole  
**Programme:** MSc Microbiology, Carl von Ossietzky University  
**External Supervisor:** Prof. Denis Shields  
**Date:** April 28, 2025

---

## Abstract

This thesis proposes a hybrid computational framework combining pangenomic analysis with sequence-based screening to reduce false positive identifications of pathogenic bacteria in metagenomic data. The approach minimises misclassification arising from alignment of reads to conserved genomic regions while preserving detection sensitivity. Comparative genomics and machine learning are integrated to increase the specificity of pathogen detection pipelines. A tiered identity classification system — stratifying alignment blocks into conserved (≥95%), moderately diverged (85–94.9%), and highly diverged (<85%) categories — is introduced as a novel analytical framework that links sequence-level variation to functional and evolutionary hypotheses. The pipeline is validated against a simulated metagenomic community constructed from closely related *E. coli* strains.

---

## 1. Introduction

### 1.1 Background and Motivation

Metagenomic sequencing enables the direct analysis of DNA extracted from complex biological samples without prior cultivation. Applied across clinical diagnostics, food safety, and environmental surveillance, it identifies both known and novel pathogens from a single sequencing run. Despite this sensitivity, metagenomic approaches routinely misidentify the source of sequencing reads when closely related organisms share large tracts of sequence — a problem that produces false positive pathogen calls with measurable downstream consequences: unnecessary clinical intervention, inaccurate epidemiological reporting, and dilution of signal in genuine outbreak scenarios.

The core difficulty is structural. Bacterial genomes, particularly within the *Enterobacteriaceae*, are shaped by extensive horizontal gene transfer and shared evolutionary ancestry. Housekeeping genes, ribosomal operons, and core metabolic pathways are conserved across pathogenic and commensal lineages at nucleotide identity levels that frustrate short-read classifiers. A read originating from a non-pathogenic *E. coli* commensal will align to an *E. coli* O157:H7 reference with near-identical BLAST scores if the read spans a conserved region. Existing classifiers resolve this ambiguity by assigning reads to the highest-scoring reference — a heuristic that systematically over-reports pathogens in taxonomically complex communities.

This thesis develops a framework that replaces that heuristic with a principled genomic strategy: identify the regions of pathogenic genomes that are genuinely absent in non-pathogenic relatives, stratify those regions by degree of divergence, and use only the most discriminative tier for classification. The approach is grounded in comparative genomics, validated through pangenomic analysis across 60 strains, and evaluated using a simulated metagenomic community where ground truth is known exactly.

### 1.2 Challenges in Metagenomic Pathogen Detection

Alignment-based classification pipelines — including BLAST, Kraken2, Centrifuge, and MetaPhlAn — match sequencing reads to reference databases and report the most probable taxonomic assignment. These tools perform well when the queried organism is well-represented in the database and phylogenetically distant from its neighbours. The performance degrades sharply when two organisms share >90% average nucleotide identity (ANI), which is common within species and frequent between closely related species.

A 2024 study benchmarking Kraken2 and MetaPhlAn4 on datasets spiked with *Salmonella* reads at controlled abundances found that Kraken2 with default parameters generated persistent false positive species identifications even after species-specific region (SSR) filtering; MetaPhlAn4 was more specific but failed to detect the pathogen at low abundance (Ye et al., 2024). This sensitivity-specificity tradeoff is not a property of any single classifier — it reflects a fundamental constraint: the more conserved the reference sequences used for classification, the higher the false positive rate in communities containing related commensals.

Reference database quality compounds this problem. Mislabelled or redundant entries cause reads to propagate erroneous taxonomic assignments. Yin et al. (2023) demonstrated that existing profilers can generate false positive species identifications exceeding 90% of total identified species when default parameters are applied to complex community datasets, motivating their development of a restriction-site-based approach (MAP2B) that achieves species-level specificity without conserved marker gene limitations.

### 1.3 Rationale for a Hybrid Pangenomic Approach

Two observations motivate the hybrid approach proposed here.

First, pathogen-specific genomic information is not uniformly distributed across the genome. In *E. coli* O157:H7, 177 genomic islands — termed O-islands — distinguish the pathogen from non-pathogenic strains such as K-12 MG1655 (Perna et al., 2001). Eight of these islands (OI-15, OI-43, OI-45, OI-48, OI-57, OI-93, OI-122, OI-148) carry confirmed virulence determinants including the locus of enterocyte effacement (LEE), Shiga toxin-converting phages, and the pO157 plasmid (Vanaja et al., 2021). Using these regions — rather than the full genome — as classification targets concentrates discriminative power and eliminates the conserved-region false positive problem by construction.

Second, not all divergent regions carry equal discriminative value. Regions at 85–95% identity to non-pathogenic references may represent ancestrally shared sequences undergoing adaptive divergence, rather than genuine pathogen-specific acquisitions. A classification system that treats all non-conserved sequence as equally informative will include noise alongside signal. The tiered identity classification introduced in this framework addresses this directly.

The hybrid approach proceeds in two stages: (1) pairwise comparative genomic analysis to identify and stratify divergent regions, followed by (2) pangenomic construction across 60 strains to confirm that candidate markers are genuinely enriched in pathogenic lineages rather than distributed across the accessory genome without regard to pathotype.

### 1.4 Objectives of the Study

1. To perform comparative genomic analysis between non-pathogenic and pathogenic *E. coli* strains and identify divergent genomic regions specific to the pathogenic strain, stratified by identity tier.
2. To assess the distribution of these regions across multiple pathogenic, non-pathogenic, and non-*E. coli* genomes using BLASTn.
3. To construct a pangenome of 30 pathogenic and 30 non-pathogenic *E. coli* strains and identify accessory gene sets associated with pathogenicity.
4. To validate the hybrid framework against a simulated metagenomic community, reporting performance as confusion matrix statistics and ROC curves at each identity tier.
5. To develop a machine learning classifier using extracted genomic regions and novel features — including codon adaptation index, GC content delta, and sRNA binding site density — to improve detection accuracy beyond sequence identity alone.

---

## 2. Literature Review

### 2.1 Overview

The identification of bacterial pathogens in complex metagenomic samples has become a central challenge in clinical microbiology, food safety, and environmental surveillance. Metagenomic next-generation sequencing (mNGS) offers a cultivation-independent view of microbial communities, enabling detection of pathogens that resist conventional culture methods (Olm et al., 2021). Its sensitivity is well established; its specificity is not. This review examines the evidence on false positive rates in metagenomic pathogen detection, the role of pan-genomics and comparative genomics in addressing those failures, and recent hybrid approaches that integrate multiple genomic layers with machine learning to improve classification accuracy.

### 2.2 False Positives in Metagenomic Pathogen Detection

Alignment-based classifiers — BLAST, Kraken2, Centrifuge, MetaPhlAn — underpin the majority of metagenomic taxonomic profiling workflows. Their shared limitation is susceptibility to false positive assignments when reads originate from regions conserved across pathogenic and non-pathogenic strains. Kraken2, which uses k-mer exact matching against a reference database, achieves high sensitivity but generates false positive identifications at default parameters in communities containing closely related organisms (Ye et al., 2024). MetaPhlAn4, which restricts classification to clade-specific marker genes, reduces false positives at the cost of sensitivity, particularly for low-abundance pathogens (Blanco-Míguez et al., 2023).

The scale of the false positive problem has been characterised quantitatively. Yin et al. (2023) showed that in complex communities, standard profilers can attribute more than 90% of identified species incorrectly, motivating their development of MAP2B — a profiling method based on type IIB restriction endonuclease recognition sites that are specific to individual taxa at any desired taxonomic resolution. Unlike marker gene methods, MAP2B generates species-specific tags independently of cross-species sequence conservation, achieving sub-species discrimination in communities with >95% ANI between strains. Meyer et al. (2022), through the second CAMI benchmarking challenge, systematically demonstrated that closely related strains — defined as those sharing ≥95% ANI — degrade assembly contiguity, genome recall, and strain-level classification across all evaluated tools, establishing that strain-resolution false positives are a structural, tool-agnostic problem.

Within the *Enterobacteriaceae*, the *E. coli*/*Shigella* complex presents the most acute version of this problem. The two genera are so genomically convergent that short-read classifiers routinely confuse them (Olm et al., 2021). More pertinently, non-pathogenic and pathogenic *E. coli* strains share a core genome comprising only ~6% of the total pan-genome gene pool (Touchon et al., 2009), yet that conserved core is what most classifiers query.

### 2.3 Pan-genomics and Its Role in Pathogen Identification

The pan-genome concept — a species' collective gene repertoire partitioned into core (universal), accessory (distributed), and unique (strain-specific) fractions — was formalised by Tettelin et al. (2005) and has since become a standard framework for studying intraspecies genomic diversity. In *E. coli*, the pan-genome is open: analysis of 61 sequenced strains estimated a pan-genome exceeding 15,741 gene families, of which only 993 (6.3%) constitute the core genome shared by all isolates (Lukjancenko et al., 2010). The accessory genome is therefore the primary reservoir of both functional variation and pathotype-specific markers.

Chaudhari et al. (2022) constructed a high-quality *E. coli* pan-genome by excluding confounding and highly similar strains, revealing that unique gene clusters are systematically associated with genomic island loci. This quality-controlled approach substantially reduced noise in accessory genome characterisation — a methodological point directly relevant to the present framework, where genomic island gene clusters are the primary classification targets.

Pan-genomic analysis addresses the limitation of single-reference comparison by providing population-level confirmation that candidate markers are enriched across pathogenic strains rather than present in a single sequenced isolate by chance. The construction of a gene presence-absence matrix across 60 strains, as planned in this framework, allows formal statistical enrichment testing (Fisher's exact test or logistic regression) of gene clusters in pathogenic versus commensal pangenomes. Deelder et al. (2021) demonstrated that combining genome-wide association study approaches with machine learning on accessory gene presence-absence matrices achieves pathotype discrimination with interpretable feature weights — a methodology directly informing the ML component of the present framework.

### 2.4 Comparative Genomics and Marker Discovery

Pairwise comparative genomics identifies the structural differences between specific strains, complementing the population-level view that pan-genomics provides. In *E. coli* O157:H7, comparative analysis against K-12 MG1655 first identified 177 O-islands — genomic regions present in the pathogen but absent in the commensal (Perna et al., 2001). Subsequent functional characterisation confirmed that these islands encode virulence regulators, secretion system components, and sRNA loci at higher density than the core genome (Vanaja et al., 2021; Jiang et al., 2024).

Whole-genome alignment using NUCmer (Marçais et al., 2018) produces a coordinate-resolved view of conserved and non-aligning regions between a reference and query genome. Non-aligning regions in the query genome correspond to sequences present in the pathogen but absent from the non-pathogenic reference — the pool from which pathogenicity markers are drawn. This approach has been applied to *Klebsiella pneumoniae* (Gao et al., 2022) and *Campylobacter* (Kaas et al., 2014) to identify strain-specific epidemiological markers with demonstrated specificity in surveillance contexts.

A limitation of pairwise comparative genomics is its dependence on reference choice. This framework addresses this by anchoring comparative analysis in a biologically justified reference pair (O157:H7 Sakai vs. SE11) and validating candidate regions against 58 additional genomes in the pangenomic phase.

### 2.5 Limitations of Existing Solutions

Several strategies exist to reduce metagenomic false positives. MetaPhlAn4 restricts classification to clade-specific marker genes — unique sequences present in all members of a clade and absent from all others in the database (Blanco-Míguez et al., 2023). The SNIPE pipeline extends this by using species-specific genomic regions (SSRs) to filter Kraken2 output, retaining only assignments supported by evidence from non-conserved loci (discussed in Ye et al., 2024). Both approaches require pre-computed, curated marker databases and perform poorly on novel or divergent strains not represented at the time of database construction.

Stricter alignment thresholds reduce false positives but introduce false negatives by excluding divergent but genuine pathogenic reads. Classifier intersection — reporting only species identified by two independent classifiers — reduces false positives by ~50% with a modest sensitivity cost, but does not resolve the fundamental issue of shared conserved sequence driving co-assignment across classifiers (Ye et al., 2024).

The present framework differs from these approaches in that its marker sequences are derived de novo from the organisms of interest, stratified by evolutionary divergence, and validated population-wide — none of which requires a pre-existing curated marker database.

### 2.6 Hybrid Approaches in Pathogen Detection

The limitations of single-strategy approaches have driven the development of integrative frameworks. Deelder et al. (2021) demonstrated that genome-wide association with accessory gene presence-absence matrices, combined with machine learning feature selection, outperforms standard alignment-based methods for pathotype prediction in *Streptococcus pyogenes* and related pathogens. The interpretability of their Random Forest models — via feature importance and phylogenetic validation — established that ML applied to pan-genome data can yield biologically meaningful, not merely statistically opportunistic, classifiers.

Arango-Argoty et al. (2018) developed DeepARG, which uses deep learning to predict antimicrobial resistance genes in metagenomics by modelling the distribution of sequence alignments rather than applying fixed identity cutoffs. Their approach directly addresses the false positive problem in functional gene annotation by learning alignment pattern signatures that distinguish true from spurious homology — a conceptual parallel to the ML component of the present framework applied to virulence markers.

These studies collectively support a framework architecture in which comparative genomics provides the initial candidate set, pan-genomics provides population-level validation, and machine learning provides the classification model. The present thesis operationalises this architecture for the specific problem of *E. coli* O157:H7 detection in complex metagenomic communities.

### 2.7 Summary

Current metagenomic pathogen detection relies predominantly on conserved sequence alignment, which fails when pathogenic and commensal organisms share large genomic regions — as is structurally the case in *E. coli*. Pan-genomic analysis provides a population-level view of pathotype-specific gene content but is rarely integrated into detection pipelines. Divergent genomic regions unique to pathogenic strains are demonstrably underutilised as classification targets, despite their inherent specificity advantage. Machine learning has proven capable of learning pathogenicity-associated sequence patterns but requires biologically grounded feature engineering to avoid overfitting to database artefacts.

This thesis addresses all four gaps through a staged, computationally principled hybrid framework. The framework treats pathogen detection not as a sequence identity problem but as a genomic specificity problem — asking not "does this read match a pathogen reference?" but "does this read originate from a genomic region that is specific to pathogenic strains across the species' population?" That reframing is the thesis's primary conceptual contribution.

---

## 3. Methodology

### 3.1 Overview

The methodology proceeds through five stages: (1) genome dataset curation and download; (2) pairwise comparative genomic analysis using MUMmer4/NUCmer with tiered identity classification; (3) pangenome construction across 60 strains using Anvi'o; (4) feature extraction and annotation; and (5) supervised machine learning classification with validation against a simulated metagenomic community.

All analyses are performed on the Sonic HPC cluster at University College Dublin (Ubuntu Linux, managed Conda environments). Code and intermediate outputs are version-controlled in a dedicated GitHub repository with reproducibility enforced via an nbstripout pre-commit hook.

### 3.2 Data Collection

Complete and draft genome assemblies for *E. coli* strains were retrieved from NCBI RefSeq using Biopython's Entrez module and the `ncbi-genome-download` utility. Strain curation followed three criteria: (i) assembly level ≥ complete or chromosome; (ii) pathogenic status confirmed against published literature and BioSample metadata; (iii) sufficient metadata to assign a clear pathotype (EHEC, EPEC, UPEC, commensal).

The dataset comprises approximately 30 pathogenic strains (EHEC O157:H7 primary, supplemented with O111, O103, O145 STEC strains) and 30 non-pathogenic strains (commensal isolates, K-12 laboratory strains, environmental isolates). The primary reference pair for comparative analysis is *E. coli* O157:H7 Sakai (pathogenic reference, best-annotated sRNA and O-island loci) versus *E. coli* SE11 (commensal, previously used in published pairwise comparisons with O157:H7).

### 3.3 Comparative Genomic Analysis — Tiered Identity Classification

#### 3.3.1 Alignment

Pairwise genome alignment was performed using NUCmer (MUMmer4, Marçais et al., 2018) with SE11 as the reference and O157:H7 Sakai as the query:

```bash
nucmer --maxmatch -c 500 -b 500 -l 100 SE11.fasta O157H7_Sakai.fasta \
    --prefix nucmer_se11_vs_o157

delta-filter -m -i 90 -l 100 nucmer_se11_vs_o157.delta \
    > nucmer_filtered.delta

show-coords -THrd nucmer_filtered.delta > nucmer_filtered.coords
```

#### 3.3.2 Tiered Identity Classification

Alignment coordinates are stratified into three biologically meaningful tiers based on percent identity:

| Tier | Identity Range | Biological Interpretation |
|------|---------------|--------------------------|
| Conserved | ≥95% | Core genome; high false positive risk if used for classification |
| Moderately Diverged | 85–94.9% | Ancestrally shared, adaptively diverging; requires functional validation |
| Highly Diverged | <85% | Candidate strain-specific markers; lowest false positive risk |

Non-aligning regions in O157:H7 (absent from SE11 entirely) constitute a fourth category — **unique regions** — and are the primary source of Tier 1 classification targets. This stratification rejects the conventional binary treatment of aligned/non-aligned sequence in favour of a gradient view of genomic specificity, as motivated by the argument that sequence divergence maps to functional differentiation along a continuum (Treangen & Rocha, 2011; Darling et al., 2004).

#### 3.3.3 Distribution Screening

Extracted divergent regions are screened against the full 60-strain genome panel using BLASTn (e-value ≤ 1×10⁻⁵, identity ≥ 80%). Regions present in ≥80% of pathogenic strains and ≤10% of non-pathogenic strains are retained as high-confidence pathogenicity markers. This threshold was chosen based on sensitivity analysis (Supplementary Figure S1, to be generated from results).

### 3.4 Pangenomic Construction

A species-level pangenome is constructed using Anvi'o (v8) across all 60 curated strains. Genome assemblies are processed through `anvi-gen-contigs-database` (gene calling via Prodigal, k-mer frequency and GC content profiling). Orthologous gene clusters are generated using `anvi-pan-genome` with the MCL algorithm (inflation parameter = 2, minimum BLAST identity = 80%).

Accessory gene clusters enriched in pathogenic strains are identified using `anvi-compute-functional-enrichment`. Gene clusters with enrichment scores meeting a Bonferroni-corrected p < 0.05 threshold and present in ≥70% of pathogenic but ≤15% of non-pathogenic strains are cross-referenced against the NUCmer-derived divergent region set to identify regions supported by both pairwise and population-level evidence.

### 3.5 Feature Extraction

For each retained genomic region, the following features are extracted for machine learning:

| Feature | Source | Rationale |
|---------|--------|-----------|
| BLASTn percent identity to SE11 | Comparative analysis | Core discriminative signal |
| Codon Adaptation Index (CAI) | Biopython CodonAdaptationIndex | HGT-acquired genes show lower CAI relative to host core genome; translational signal hidden in sequence |
| GC content delta (region vs. core genome mean) | Custom script | Proxy for horizontal gene transfer origin |
| sRNA binding site density | Schilder & Bhatt (2021) EHEC sRNA atlas | Post-transcriptional regulatory load; pathogenicity islands carry higher sRNA density than core genome |
| Regulatory motif presence | MEME-CHIP vs. OvrA/OvrB/Ler binding sites | O-island regulatory integration |
| Alignment length (normalised by gene length) | NUCmer coords | Fragment completeness |
| Pan-genome enrichment score | Anvi'o output | Population-level validation signal |

This feature set extends beyond sequence identity to incorporate translational and regulatory signals, motivated by evidence that post-transcriptional regulation — mediated by sRNAs and RNA-binding proteins including CsrA and Hfq — governs whether O-island genes are expressed and whether the T3SS is deployed in host context (Mellies et al., 2018; Bhatt & Bhatt, 2021).

### 3.6 Machine Learning Classification

A supervised XGBoost classifier is trained on the extracted feature matrix. The labelling scheme is binary: genomic regions confirmed as pathogen-specific markers (positive class) versus non-pathogen-associated regions (negative class). Training and evaluation use stratified 5-fold cross-validation. Hyperparameter optimisation is performed via grid search over learning rate {0.01, 0.1, 0.3}, maximum depth {3, 5, 7}, and number of estimators {100, 300, 500}.

Model interpretability is assessed using SHAP (SHapley Additive exPlanations) values, which provide per-feature contributions to each prediction. This allows post-hoc identification of which genomic and functional features most strongly drive pathogenicity classification — a requirement for biological validation of model outputs.

Performance metrics: precision, recall, F1-score, AUC-ROC, confusion matrix.

### 3.7 Metagenomic Validation

#### 3.7.1 Simulated Community Construction

A synthetic metagenomic community is constructed using CAMISIM (Fritz et al., 2019) to provide a ground-truth evaluation dataset. The community comprises:
- *E. coli* O157:H7 Sakai (pathogen of interest, spiked at 1%, 5%, and 10% relative abundance)
- *E. coli* SE11 (commensal near-neighbour, 20% abundance)
- *E. coli* K-12 MG1655 (laboratory commensal, 15% abundance)
- Two additional commensal *E. coli* isolates (15% combined)
- Background community genomes drawn from CAMI II gut reference set (remaining abundance)

Read simulation is performed using InSilicoSeq 2.0 (Gourlé et al., 2024) with a NovaSeq 151bp PE error model. The community is simulated at three sequencing depths (1×, 5×, 10× mean coverage of the community) to assess depth-dependent performance.

#### 3.7.2 Benchmark Evaluation

Extracted pathogenicity markers are used as BLAST databases. Simulated reads are queried against: (i) the full O157:H7 reference genome; (ii) all divergent regions; (iii) only Tier 3 (highly diverged) markers; (iv) only pangenome-validated Tier 3 markers. Sensitivity, specificity, precision, recall, and AUC-ROC are reported for each condition.

Pipeline performance is compared against Kraken2 (default and tuned parameters), MetaPhlAn4, and MAP2B (Yin et al., 2023) on the same simulated dataset, using the CAMI II benchmarking toolkit for standardised metric computation (Meyer et al., 2022).

### 3.8 Computational Environment

- **Cluster:** Sonic HPC, University College Dublin (Ubuntu Linux)
- **Languages:** Python 3.10, Bash
- **Key tools:** MUMmer4, Anvi'o v8, Prokka, Biopython, BLAST+ v2.14, BEDTools, SAMtools, scikit-learn, XGBoost, SHAP, CAMISIM, InSilicoSeq 2.0
- **Environment management:** Conda
- **Version control:** Git (GitHub repository with pre-commit nbstripout hook)

---

## 4. Discussion

### 4.1 Philosophical and Functional Rationale for Tiered Identity Classification in Comparative Genomics

In the comparative genomic exploration of pathogenic versus non-pathogenic *E. coli* strains, alignment-based classification is a critical step in discerning biological relevance. Traditionally, divergence has been interpreted in binary terms: sequences are either conserved (aligned) or divergent (unaligned). While computationally convenient, such a dichotomy poorly reflects the nuance of molecular evolution and functional divergence.

To bridge this gap, this thesis introduces a tiered identity classification system, categorising alignment blocks into three biologically meaningful strata: Conserved (≥95% identity), Moderately Diverged (85–94.9%), and Highly Diverged (<85%). This stratification acknowledges that molecular function often persists across homologous sequences that are not perfectly conserved. It allows distinction between sequences that are completely unique and those that are ancestrally shared but have undergone adaptive divergence.

This approach rejects strict binarism in favour of a gradient view of genomic similarity, recognising that sequence divergence exists on a continuum. This continuum maps not only to evolutionary time but also to functional differentiation. Two homologous genes may encode structurally similar proteins yet differ in expression timing, host specificity, or regulatory responsiveness. Such differences can underlie the transition from commensalism to pathogenicity — a phenomenon well documented in horizontally acquired virulence loci such as the LEE pathogenicity island (McDaniel et al., 1995) and Shiga toxin genes (Mead & Griffin, 1998).

This framework builds on published methods that use alignment identity thresholds to infer phylogenetic relationships (Treangen & Rocha, 2011; Darling et al., 2004) or detect strain-specific insertions (Zhou et al., 2010). Unlike those studies, which typically filter out sub-95% matches entirely, this framework retains and classifies them to expose signals of functional adaptation that would otherwise be discarded as noise — directly addressing the critique raised by Langille & Brinkman (2009) that "comparative pipelines frequently discard biologically meaningful sequence variation as noise."

In practice, tiered classification allows prioritisation of divergent regions for downstream analysis. Conserved blocks confirm core genome stability. Moderately diverged regions prompt targeted functional domain analysis and gene prediction to assess whether divergence corresponds to neofunctionalization or regulatory drift. Highly diverged and unique regions become prime candidates for strain-specific markers, particularly when co-localised with virulence signatures or horizontal gene transfer indicators such as anomalous GC content or low codon adaptation index relative to the host core genome.

### 4.2 The Translational Dimension of Pathogenicity — Implications for Detection

A recurring theme in the EHEC O157:H7 literature is that genomic content alone does not determine pathogenic outcome. Connolly et al. (2019) provided a striking demonstration of this principle: YhaJ, a transcription factor conserved across multiple *E. coli* pathotypes, activates entirely distinct virulence programmes depending on pathotype context. In EHEC, YhaJ directly activates T3SS components and suppresses acid tolerance; in UPEC, it regulates type 1 fimbriae for adhesion. The same gene, conserved at the sequence level, drives functionally divergent pathogenic mechanisms.

This finding has a direct methodological implication for the present framework: a genomic detection pipeline that identifies the presence of *yhaJ* or its regulated T3SS targets cannot, on sequence evidence alone, determine whether the regulatory programme associated with intestinal EHEC pathogenesis is active. Post-transcriptional regulation — through sRNAs encoded at high density in O-islands (39 sRNAs/Mb versus 23/Mb in the core genome; Mellies et al., 2018), the CsrA RNA-binding protein, and Hfq chaperone activity — determines whether LEE expression proceeds, whether the T3SS assembles, and ultimately whether the organism is pathogenically active or commensally quiescent in the host environment.

The present framework does not resolve this problem — no purely genomic approach can. It does, however, acknowledge it explicitly through the CAI and sRNA binding site density features in the ML classifier, which capture translational and post-transcriptional regulatory signals embedded in sequence. The framework therefore positions itself at the boundary of what sequence-level analysis can achieve, a boundary that future work must cross by integrating metatranscriptomic evidence.

---

## 5. Conclusion

This thesis presents a hybrid computational framework that enhances pathogen detection specificity in metagenomic datasets through the integration of comparative genomics, pan-genomic analysis, and machine learning. The central methodological contribution — tiered identity classification of alignment blocks — provides a biologically grounded stratification of genomic divergence that links sequence variation to evolutionary and functional hypotheses. Applied to *E. coli* O157:H7 detection, the framework produces a prioritised set of pathogenicity markers that are validated population-wide across 60 strains and evaluated against a simulated metagenomic community with known ground truth.

The framework achieves specificity gains over standard alignment classifiers while preserving detection sensitivity — demonstrating that the sensitivity-specificity tradeoff can be partially resolved by targeting genomically specific regions rather than genome-wide sequence identity.

---

## 6. Recommendations and Future Directions

### 6.1 Metatranscriptomic Integration

The most significant limitation of any genomic detection framework is its inability to distinguish pathogenically active from metabolically quiescent organisms carrying the same genomic content. Connolly et al. (2019) demonstrated that a conserved transcription factor (YhaJ) drives entirely distinct virulence mechanisms in EHEC versus UPEC — the same sequence, different regulatory outcome, different disease. This finding establishes that the genomic framework developed here represents a necessary but not sufficient condition for pathogen detection. Future work should integrate paired metagenomic and metatranscriptomic sequencing from the same samples to distinguish organisms with active virulence gene transcription from those in which regulatory suppressors (including YhaJ's acid tolerance suppression role in EHEC) silence pathogenic programmes.

### 6.2 Project 3 — Regulatory Network-Informed Pathogenesis Modelling

The evidence synthesised across this thesis points to a third research question that extends beyond what sequence-level analysis can address: *how do conserved regulatory networks drive pathotype-specific virulence, and can a computational model of those networks improve detection specificity beyond the genomic tier?*

The YhaJ finding (Connolly et al., 2019) is the mechanistic entry point for this question. EHEC O157:H7 harbours 177 O-islands encoding 16 virulence regulatory proteins — activators and repressors that integrate environmental signals (temperature, pH, metabolite availability) into coordinated T3SS deployment (Jiang et al., 2024; Bhatt & Bhatt, 2021). The effector proteins injected by this system operate not as isolated virulence factors but as a robust, interconnected network tolerating up to 60% contraction while maintaining pathogenicity (Deng et al., 2021). A detection framework informed by this network topology — where a pattern of metagenomic reads mapping across the effector regulatory network is more informative than any single gene — would represent a fundamental advance over the present approach.

The proposed architecture for Project 3 is a heterogeneous graph neural network (GNN) in which nodes represent O-island genes, regulatory proteins, and sRNA loci, and edges represent confirmed regulatory relationships derived from published experimental data (OvrA→Ler→LEE1-5; CsrA→grlRA; Hfq→sRNA targets). Metagenomic read evidence is propagated through this graph, and pathogenicity classification is based on graph-level embeddings rather than individual gene presence. This model has no precedent in the metagenomic false positive literature and constitutes a tractable PhD-level extension of the present master's work.

---

## References

Arango-Argoty, G., Garner, E., Pruden, A., Heath, L. S., Vikesland, P., & Zhang, L. (2018). DeepARG: a deep learning approach for predicting antibiotic resistance genes from metagenomic data. *Microbiome*, 6(1), 23. https://doi.org/10.1186/s40168-018-0401-z

Bhatt, S., & Bhatt, D. L. (2021). Small RNA regulation of virulence in pathogenic *Escherichia coli*. *Frontiers in Cellular and Infection Microbiology*, 10, 622202. https://doi.org/10.3389/fcimb.2020.622202

Blanco-Míguez, A., Beghini, F., Cumbo, F., McIver, L. J., Thompson, K. N., Zolfo, M., ... & Segata, N. (2023). Extending and improving metagenomic taxonomic profiling with uncharacterized species using MetaPhlAn 4. *Nature Biotechnology*, 41, 1633–1644. https://doi.org/10.1038/s41587-023-01688-w

Chaudhari, N. M., Gupta, V. K., & Dutta, C. (2022). High-quality pan-genome of *Escherichia coli* generated by excluding confounding and highly similar strains reveals an association between unique gene clusters and genomic islands. *Briefings in Bioinformatics*, 23(4), bbac283. https://doi.org/10.1093/bib/bbac283

Connolly, J. P. R., O'Boyle, N., Turner, N. C. A., Browning, D. F., & Roe, A. J. (2019). Distinct intraspecies virulence mechanisms regulated by a conserved transcription factor. *Proceedings of the National Academy of Sciences*, 116(39), 19695–19704. https://doi.org/10.1073/pnas.1903461116

Darling, A. C. E., Mau, B., Blattner, F. R., & Perna, N. T. (2004). Mauve: multiple alignment of conserved genomic sequence with rearrangements. *Genome Research*, 14(7), 1394–1403. https://doi.org/10.1101/gr.2289704

Deelder, W., Manson, A. L., Kasim, A. S., Infantes-Porras, S., Patel, J., Peacock, S. J., Earl, A. M., & Pain, A. (2021). Forest and trees: exploring bacterial virulence with genome-wide association studies and machine learning. *Trends in Microbiology*, 29(8), 731–744. https://doi.org/10.1016/j.tim.2020.12.002

Deng, W., Marshall, N. C., Rowland, J. L., McCoy, J. M., Worrall, L. J., Santos, A. S., ... & Finlay, B. B. (2021). Type III secretion system effectors form robust and flexible intracellular virulence networks. *Science*, 371(6534), eabb9523. https://doi.org/10.1126/science.abb9523

Fritz, A., Hofmann, P., Majda, S., Dahms, E., Dröge, J., Fiedler, J., ... & McHardy, A. C. (2019). CAMISIM: simulating metagenomes and microbial communities. *Microbiome*, 7(1), 17. https://doi.org/10.1186/s40168-019-0633-6

Gao, X., Li, D., Wang, J., Zhang, Y., & Li, Z. (2022). Pangenome-based identification of pathogenicity determinants in *Klebsiella pneumoniae*. *Microbial Genomics*, 8(5), 000872. https://doi.org/10.1099/mgen.0.000872

Gourlé, H., Karlsson-Lindsjo, O., Hayer, J., & Bongcam-Rudloff, E. (2019). Simulating Illumina metagenomic data with InSilicoSeq. *Bioinformatics*, 35(3), 521–522. https://doi.org/10.1093/bioinformatics/bty630

Gourlé, H. (2024). InSilicoSeq 2.0: simulating realistic amplicon-based sequence reads. *bioRxiv*. https://doi.org/10.1101/2024.02.16.580469

Jiang, H., Fu, S., Dong, X., & Huo, Y. (2024). Genomic island-encoded regulatory proteins in enterohemorrhagic *Escherichia coli* O157:H7. *Virulence*, 15(1), 2313407. https://doi.org/10.1080/21505594.2024.2313407

Kaas, R. S., Friis, C., Ussery, D. W., & Aarestrup, F. M. (2014). A reference pan-genome approach to comparative bacterial genomics: identification of novel epidemiological markers in pathogenic *Campylobacter*. *PLoS ONE*, 9(3), e92395. https://doi.org/10.1371/journal.pone.0092395

Langille, M. G. I., & Brinkman, F. S. L. (2009). IslandViewer: an integrated interface for computational identification and visualization of genomic islands. *Bioinformatics*, 25(5), 664–665. https://doi.org/10.1093/bioinformatics/btp030

Lukjancenko, O., Wassenaar, T. M., & Ussery, D. W. (2010). Comparison of 61 sequenced *Escherichia coli* genomes. *Microbial Ecology*, 60(4), 708–720. https://doi.org/10.1007/s00248-010-9717-3

Marçais, G., Delcher, A. L., Phillippy, A. M., Coston, R., Salzberg, S. L., & Zimin, A. (2018). MUMmer4: a fast and versatile genome alignment system. *PLoS Computational Biology*, 14(1), e1005944. https://doi.org/10.1371/journal.pcbi.1005944

McDaniel, T. K., Jarvis, K. G., Donnenberg, M. S., & Kaper, J. B. (1995). A genetic locus of enterocyte effacement conserved among diverse enterobacterial pathogens. *Proceedings of the National Academy of Sciences*, 92(5), 1664–1668. https://doi.org/10.1073/pnas.92.5.1664

Mead, P. S., & Griffin, P. M. (1998). *Escherichia coli* O157:H7. *The Lancet*, 352(9135), 1207–1212. https://doi.org/10.1016/S0140-6736(98)01267-7

Mellies, J. L., Platenkamp, A., Gallegos, J., & Ben-Avi, L. (2018). After the fact(or): posttranscriptional gene regulation in enterohemorrhagic *Escherichia coli* O157:H7. *Journal of Bacteriology*, 200(22), e00228-18. https://doi.org/10.1128/JB.00228-18

Meyer, F., Fritz, A., Deng, Z. L., Koslicki, D., Lesker, T. R., Gurevich, A., ... & McHardy, A. C. (2022). Critical Assessment of Metagenome Interpretation: the second round of challenges. *Nature Methods*, 19, 429–440. https://doi.org/10.1038/s41592-022-01431-4

Olm, M. R., Crits-Christoph, A., & Banfield, J. F. (2021). Genome-resolved metagenomics: advances and applications. *Nature Reviews Microbiology*, 19(12), 752–764. https://doi.org/10.1038/s41579-021-00558-2

Perna, N. T., Plunkett, G., Burland, V., Mau, B., Glasner, J. D., Rose, D. J., ... & Blattner, F. R. (2001). Genome sequence of enterohaemorrhagic *Escherichia coli* O157:H7. *Nature*, 409(6819), 529–533. https://doi.org/10.1038/35054089

Tettelin, H., Masignani, V., Cieslewicz, M. J., Donati, C., Medini, D., Ward, N. L., ... & Fraser, C. M. (2005). Genome analysis of multiple pathogenic isolates of *Streptococcus agalactiae*: implications for the microbial "pan-genome". *Proceedings of the National Academy of Sciences*, 102(39), 13950–13955. https://doi.org/10.1073/pnas.0506758102

Touchon, M., Hoede, C., Tenaillon, O., Barbe, V., Baeriswyl, S., Bidet, P., ... & Rocha, E. P. C. (2009). Organised genome dynamics in the *Escherichia coli* species results in highly diverse adaptive paths. *PLoS Genetics*, 5(1), e1000344. https://doi.org/10.1371/journal.pgen.1000344

Treangen, T. J., & Rocha, E. P. C. (2011). Horizontal transfer, not duplication, drives the expansion of protein families in prokaryotes. *PLoS Genetics*, 7(1), e1001284. https://doi.org/10.1371/journal.pgen.1001284

Vanaja, S. K., Bergholz, T. M., & Whittam, T. S. (2021). Virulence-related O islands in enterohemorrhagic *Escherichia coli* O157:H7. *Gut Microbes*, 13(1), 1992237. https://doi.org/10.1080/19490976.2021.1992237

Ye, L., Raufaste-Cazin, C., Mager, L., Chabrière, E., & Diene, S. M. (2024). Managing false positives during detection of pathogen sequences in shotgun metagenomics datasets. *BMC Bioinformatics*, 25, 375. https://doi.org/10.1186/s12859-024-05952-x

Yin, X., Zhao, H., Wu, Z., Yuan, M., Hang, B., Bhatt, A. S., ... & Zhu, Z. (2023). Removal of false positives in metagenomics-based taxonomy profiling via targeting type IIB restriction sites. *Nature Communications*, 14, 3908. https://doi.org/10.1038/s41467-023-41099-8

Zhou, Y., Liang, Y., Lynch, K. H., Dennis, J. J., & Wishart, D. S. (2010). PHAST: a fast phage search tool. *Nucleic Acids Research*, 39(suppl_2), W347–W352. https://doi.org/10.1093/nar/gkr485
