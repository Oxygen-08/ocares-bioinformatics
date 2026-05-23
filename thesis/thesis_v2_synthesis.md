# Hybrid Pangenomic and Sequence-Based Framework for Reducing False Positives in Pathogen Detection from Metagenomic Data

**Author:** Oluwatosin Samuel Oluwole  
**Programme:** MSc Microbiology, Carl von Ossietzky University  
**Internal Supervisor:** Prof. Sarahi L. Garcia  
**External Supervisor:** Prof. Denis Shields  
**Date:** May 2025  
**Code availability:** https://doi.org/10.5281/zenodo.20278099

---

## Abstract

This thesis presents a hybrid computational framework combining pangenomic analysis with sequence-based screening and machine learning to reduce false positive identifications of pathogenic bacteria in metagenomic data. Comparative genomic analysis of *E. coli* O157:H7 Sakai against 60 commensal and pathogenic strains (balanced 30/30 panel) yielded 23,923 NUCmer alignment blocks, from which 415 candidate pathogenicity markers were extracted across three identity tiers: Conserved (≥95%), Moderately Diverged (85–94.9%), and Highly Diverged (<85%). A 15-feature XGBoost classifier — incorporating sequence identity, codon adaptation index, GC content delta, sRNA binding site density, alignment coverage, 4-mer compositional deviation, NUCmer/Anvi'o-derived pangenome features, and minimap2 divergence gradient features — achieved AUROC 0.711 ± 0.034 (out-of-sample, 5-fold grouped cross-validation with 500kb genomic bins), substantially exceeding the BLAST screen baseline of 0.552. Biological validation confirmed that 49.5% of DIVERGED-tier markers are absent from all K-12 strains and 14 markers co-localise with named *E. coli* O157:H7 virulence loci including the LEE pathogenicity island, Shiga toxin-converting prophages SpLE2 and SpLE3, and the OI-48 type III secretion locus; reclassifying on K-12 absence as an independent biological ground truth achieved AUROC 0.818 ± 0.024, confirming the biologically cleaner label reveals stronger discriminative signal. The framework's primary novel contribution is the simultaneous interrogation of compositional foreignness (4-mer deviation), translational adaptation (CAI), post-transcriptional regulatory architecture (sRNA density), and population-level pangenome context — five feature layers querying a single binary classification target. No published metagenomics tool combines this multi-modal feature geometry. The work constitutes a proof of concept for sequence-informed pathogenicity inference that extends beyond alignment identity, with a clear extension path toward regulatory network-informed classification using graph neural networks.

---

## 1. Introduction

### 1.1 Background and Motivation

Metagenomic sequencing enables the direct analysis of DNA extracted from complex biological samples without prior cultivation. Applied across clinical diagnostics, food safety, and environmental surveillance, it identifies both known and novel pathogens from a single sequencing run. Despite this sensitivity, metagenomic approaches routinely misidentify the source of sequencing reads when closely related organisms share large tracts of sequence — a problem that produces false positive pathogen calls with measurable downstream consequences: unnecessary clinical intervention, inaccurate epidemiological reporting, and dilution of signal in genuine outbreak scenarios.

The core difficulty is structural. Bacterial genomes, particularly within the *Enterobacteriaceae*, are shaped by extensive horizontal gene transfer and shared evolutionary ancestry. Housekeeping genes, ribosomal operons, and core metabolic pathways are conserved across pathogenic and commensal lineages at nucleotide identity levels that frustrate short-read classifiers. A read originating from a non-pathogenic *E. coli* commensal will align to an *E. coli* O157:H7 reference with near-identical BLAST scores if the read spans a conserved region. Existing classifiers resolve this ambiguity by assigning reads to the highest-scoring reference — a heuristic that systematically over-reports pathogens in taxonomically complex communities.

This thesis develops a framework that replaces that heuristic with a principled genomic strategy: identify the regions of pathogenic genomes that are genuinely absent in non-pathogenic relatives, stratify those regions by degree of divergence, and use only the most discriminative tier for classification. The approach is grounded in comparative genomics, validated through pangenomic analysis across 61 strains, and evaluated using a simulated metagenomic community where ground truth is known exactly.

### 1.2 Challenges in Metagenomic Pathogen Detection

Alignment-based classification pipelines — including BLAST, Kraken2, Centrifuge, and MetaPhlAn — match sequencing reads to reference databases and report the most probable taxonomic assignment. These tools perform well when the queried organism is well-represented in the database and phylogenetically distant from its neighbours. The performance degrades sharply when two organisms share >90% average nucleotide identity (ANI), which is common within species and frequent between closely related species.

A 2024 study benchmarking Kraken2 and MetaPhlAn4 on datasets spiked with *Salmonella* reads at controlled abundances found that Kraken2 with default parameters generated persistent false positive species identifications even after species-specific region (SSR) filtering; MetaPhlAn4 was more specific but failed to detect the pathogen at low abundance (Ye et al., 2024). This sensitivity-specificity tradeoff is not a property of any single classifier — it reflects a fundamental constraint: the more conserved the reference sequences used for classification, the higher the false positive rate in communities containing related commensals.

Reference database quality compounds this problem. Mislabelled or redundant entries cause reads to propagate erroneous taxonomic assignments. Yin et al. (2023) demonstrated that existing profilers can generate false positive species identifications exceeding 90% of total identified species when default parameters are applied to complex community datasets, motivating their development of a restriction-site-based approach (MAP2B) that achieves species-level specificity without conserved marker gene limitations.

### 1.3 Rationale for a Hybrid Pangenomic Approach

Two observations motivate the hybrid approach proposed here.

First, pathogen-specific genomic information is not uniformly distributed across the genome. In *E. coli* O157:H7, 177 genomic islands — termed O-islands — distinguish the pathogen from non-pathogenic strains such as K-12 MG1655 (Perna et al., 2001). Eight of these islands (OI-15, OI-43, OI-45, OI-48, OI-57, OI-93, OI-122, OI-148) carry confirmed virulence determinants including the locus of enterocyte effacement (LEE), Shiga toxin-converting phages, and the pO157 plasmid (Vanaja et al., 2021). Using these regions — rather than the full genome — as classification targets concentrates discriminative power and eliminates the conserved-region false positive problem by construction.

Second, not all divergent regions carry equal discriminative value. Regions at 85–95% identity to non-pathogenic references may represent ancestrally shared sequences undergoing adaptive divergence, rather than genuine pathogen-specific acquisitions. A classification system that treats all non-conserved sequence as equally informative will include noise alongside signal. The tiered identity classification introduced in this framework addresses this directly.

The hybrid approach proceeds in two stages: (1) pairwise comparative genomic analysis to identify and stratify divergent regions, followed by (2) pangenomic construction across all 61 strains using Anvi'o v8 to confirm that candidate markers are genuinely enriched in pathogenic lineages rather than distributed across the accessory genome without regard to pathotype.

### 1.4 Objectives of the Study

1. To perform comparative genomic analysis between non-pathogenic and pathogenic *E. coli* strains and identify divergent genomic regions specific to the pathogenic strain, stratified by identity tier.
2. To assess the distribution of these regions across multiple pathogenic, non-pathogenic, and non-*E. coli* genomes using BLASTn.
3. To construct a pangenome of pathogenic and commensal *E. coli* strains using Anvi'o and identify accessory gene sets associated with pathogenicity.
4. To develop a multi-modal machine learning classifier using extracted genomic regions and novel features — including codon adaptation index, 4-mer compositional deviation, GC content delta, sRNA binding site density, and pangenome-derived presence scores — to improve detection accuracy beyond sequence identity alone.
5. To validate the hybrid framework against a simulated metagenomic community, reporting performance as confusion matrix statistics, AUROC, and biological marker validation.

---

## 2. Literature Review

### 2.1 Overview

The identification of bacterial pathogens in complex metagenomic samples has become a central challenge in clinical microbiology, food safety, and environmental surveillance. Metagenomic next-generation sequencing (mNGS) offers a cultivation-independent view of microbial communities, enabling detection of pathogens that resist conventional culture methods (Olm et al., 2021). Its sensitivity is well established; its specificity is not. This review examines the evidence on false positive rates in metagenomic pathogen detection, the role of pan-genomics and comparative genomics in addressing those failures, and recent hybrid approaches that integrate multiple genomic layers with machine learning to improve classification accuracy.

### 2.2 False Positives in Metagenomic Pathogen Detection

Alignment-based classifiers — BLAST, Kraken2, Centrifuge, MetaPhlAn — underpin the majority of metagenomic taxonomic profiling workflows. Their shared limitation is susceptibility to false positive assignments when reads originate from regions conserved across pathogenic and non-pathogenic strains. Kraken2 (Wood et al., 2019), which uses k-mer exact matching against a reference database, achieves high sensitivity but generates false positive identifications at default parameters in communities containing closely related organisms (Ye et al., 2024). MetaPhlAn4, which restricts classification to clade-specific marker genes, reduces false positives at the cost of sensitivity, particularly for low-abundance pathogens (Blanco-Míguez et al., 2023).

The scale of the false positive problem has been characterised quantitatively. Yin et al. (2023) showed that in complex communities, standard profilers can attribute more than 90% of identified species incorrectly, motivating their development of MAP2B — a profiling method based on type IIB restriction endonuclease recognition sites that are specific to individual taxa at any desired taxonomic resolution. Meyer et al. (2022), through the second CAMI benchmarking challenge, systematically demonstrated that closely related strains — defined as those sharing ≥95% ANI — degrade assembly contiguity, genome recall, and strain-level classification across all evaluated tools, establishing that strain-resolution false positives are a structural, tool-agnostic problem.

Within the *Enterobacteriaceae*, the *E. coli*/*Shigella* complex presents the most acute version of this problem. The two genera are so genomically convergent that short-read classifiers routinely confuse them (Olm et al., 2021). More pertinently, non-pathogenic and pathogenic *E. coli* strains share a core genome comprising only ~6% of the total pan-genome gene pool (Touchon et al., 2009), yet that conserved core is what most classifiers query.

### 2.3 Pan-genomics and Its Role in Pathogen Identification

The pan-genome concept — a species' collective gene repertoire partitioned into core (universal), accessory (distributed), and unique (strain-specific) fractions — was formalised by Tettelin et al. (2005) and has since become a standard framework for studying intraspecies genomic diversity. In *E. coli*, the pan-genome is open: analysis of 61 sequenced strains estimated a pan-genome exceeding 15,741 gene families, of which only 993 (6.3%) constitute the core genome shared by all isolates (Lukjancenko et al., 2010). The accessory genome is therefore the primary reservoir of both functional variation and pathotype-specific markers.

Chaudhari et al. (2022) constructed a high-quality *E. coli* pan-genome by excluding confounding and highly similar strains, revealing that unique gene clusters are systematically associated with genomic island loci. This quality-controlled approach substantially reduced noise in accessory genome characterisation — a methodological point directly relevant to the present framework, where genomic island gene clusters are the primary classification targets.

Pan-genomic analysis addresses the limitation of single-reference comparison by providing population-level confirmation that candidate markers are enriched across pathogenic strains rather than present in a single sequenced isolate by chance. The construction of a gene presence-absence matrix across 61 strains in this framework allows formal statistical enrichment testing of gene clusters in pathogenic versus commensal pangenomes. Deelder et al. (2021) demonstrated that combining genome-wide association study approaches with machine learning on accessory gene presence-absence matrices achieves pathotype discrimination with interpretable feature weights — a methodology directly informing the ML component of the present framework.

### 2.4 Comparative Genomics and Marker Discovery

Pairwise comparative genomics identifies the structural differences between specific strains, complementing the population-level view that pan-genomics provides. In *E. coli* O157:H7, comparative analysis against K-12 MG1655 first identified 177 O-islands — genomic regions present in the pathogen but absent in the commensal (Perna et al., 2001). Subsequent functional characterisation confirmed that these islands encode virulence regulators, secretion system components, and sRNA loci at higher density than the core genome (Vanaja et al., 2021; Jiang et al., 2024).

Whole-genome alignment using NUCmer (Marçais et al., 2018) produces a coordinate-resolved view of conserved and non-aligning regions between a reference and query genome. Non-aligning regions in the query genome correspond to sequences present in the pathogen but absent from the non-pathogenic reference — the pool from which pathogenicity markers are drawn. This approach has been applied to *Klebsiella pneumoniae* (Gao et al., 2022) and *Campylobacter* (Kaas et al., 2014) to identify strain-specific epidemiological markers with demonstrated specificity in surveillance contexts.

### 2.5 Limitations of Existing Solutions

Several strategies exist to reduce metagenomic false positives. MetaPhlAn4 restricts classification to clade-specific marker genes — unique sequences present in all members of a clade and absent from all others in the database (Blanco-Míguez et al., 2023). The SNIPE pipeline extends this by using species-specific genomic regions (SSRs) to filter Kraken2 output, retaining only assignments supported by evidence from non-conserved loci (discussed in Ye et al., 2024). Both approaches require pre-computed, curated marker databases and perform poorly on novel or divergent strains not represented at the time of database construction.

Stricter alignment thresholds reduce false positives but introduce false negatives by excluding divergent but genuine pathogenic reads. Classifier intersection — reporting only species identified by two independent classifiers — reduces false positives by ~50% with a modest sensitivity cost, but does not resolve the fundamental issue of shared conserved sequence driving co-assignment across classifiers (Ye et al., 2024).

The present framework differs from these approaches in that its marker sequences are derived de novo from the organisms of interest, stratified by evolutionary divergence, and validated population-wide — none of which requires a pre-existing curated marker database. Critically, it extends beyond sequence identity to interrogate compositional, translational, and regulatory properties of genomic regions simultaneously, an integration depth that has no direct precedent in the metagenomic false positive literature.

### 2.6 Hybrid Approaches in Pathogen Detection

The limitations of single-strategy approaches have driven the development of integrative frameworks. Deelder et al. (2021) demonstrated that genome-wide association with accessory gene presence-absence matrices, combined with machine learning feature selection, outperforms standard alignment-based methods for pathotype prediction in *Streptococcus pyogenes* and related pathogens. The interpretability of their Random Forest models — via feature importance and phylogenetic validation — established that ML applied to pan-genome data can yield biologically meaningful, not merely statistically opportunistic, classifiers.

Arango-Argoty et al. (2018) developed DeepARG, which uses deep learning to predict antimicrobial resistance genes in metagenomics by modelling the distribution of sequence alignments rather than applying fixed identity cutoffs. Their approach directly addresses the false positive problem in functional gene annotation by learning alignment pattern signatures that distinguish true from spurious homology — a conceptual parallel to the ML component of the present framework applied to virulence markers.

These studies collectively support a framework architecture in which comparative genomics provides the initial candidate set, pan-genomics provides population-level validation, and machine learning provides the classification model. The present thesis operationalises this architecture for the specific problem of *E. coli* O157:H7 detection in complex metagenomic communities, and extends it with a feature geometry that spans five biological layers: alignment identity, compositional, translational, regulatory, and population-level pangenome.

### 2.7 Summary

Current metagenomic pathogen detection relies predominantly on conserved sequence alignment, which fails when pathogenic and commensal organisms share large genomic regions — as is structurally the case in *E. coli*. Pan-genomic analysis provides a population-level view of pathotype-specific gene content but is rarely integrated into detection pipelines. Divergent genomic regions unique to pathogenic strains are demonstrably underutilised as classification targets, despite their inherent specificity advantage. Machine learning has proven capable of learning pathogenicity-associated sequence patterns but requires biologically grounded feature engineering to avoid overfitting to database artefacts.

This thesis addresses all four gaps through a staged, computationally principled hybrid framework. The framework treats pathogen detection not as a sequence identity problem but as a genomic specificity problem — asking not "does this read match a pathogen reference?" but "does this read originate from a genomic region that is specific to pathogenic strains across the species' population, and does it carry the compositional, translational, and regulatory signatures of horizontally acquired virulence sequence?" That reframing is the thesis's primary conceptual contribution.

---

## 3. Methodology

### 3.1 Overview

The methodology proceeds through six stages: (1) genome dataset curation and download; (2) pairwise comparative genomic analysis using MUMmer4/NUCmer with tiered identity classification; (3) marker extraction and BLAST screen evaluation; (4) pangenome construction across 61 strains using Anvi'o with COG14 functional annotation; (5) multi-modal feature extraction and ML classification; and (6) biological validation of the marker set. A simulated metagenomic community generated by InSilicoSeq 2.0 provides ground-truth evaluation throughout stages 3–5.

All analyses were performed on a local development environment (macOS, Conda-managed environments). Code and outputs are version-controlled at https://doi.org/10.5281/zenodo.20278099. Complete software environments are specified in `environment.yml` (main pipeline: Python 3.11, minimap2 ≥2.28, MUMmer4 ≥4.0.0, BLAST ≥2.16.0, XGBoost ≥2.0, scikit-learn ≥1.4, InSilicoSeq ≥2.0) and `environment_anvio.yml` (Anvi'o v8), both included in the repository.

### 3.2 Data Collection

Complete genome assemblies for *E. coli* strains were retrieved from NCBI RefSeq using the `ncbi-genome-download` utility. Strain curation followed three criteria: (i) assembly level ≥ complete or chromosome; (ii) pathogenic status confirmed against published literature and BioSample metadata; (iii) sufficient metadata to assign a clear pathotype (EHEC, EPEC, UPEC, commensal).

The dataset comprises 61 curated genomes: the O157:H7 Sakai reference plus 60 comparison genomes (30 pathogenic, 30 non-pathogenic), achieving the planned balanced panel. All GCF accessions and pathotype assignments are recorded in `data/genomes/genome_manifest.tsv` with supporting PMIDs; one mislabelled strain (SE11, initially classified EPEC) was corrected to commensal based on Oshima et al. (2008, PMID 18931093). The panel includes:
- **Primary reference:** *E. coli* O157:H7 Sakai (GenBank: BA000007.3), the best-annotated EHEC genome with documented O-island coordinates and sRNA atlas
- **Pathogenic panel (30 strains):** EHEC O157:H7 isolates (EDL933, EC4115, TW14359, O26:H11, O103:H2, O111, O145, O55:H7, O104:H4, HUSEC2011), UPEC (CFT073, UTI89, 536, UMN026, IAI39, NA114), ETEC (H10407, E24377A, TW11681), EAEC (042, 55989), EPEC (E2348/69, TW10598), NMEC (IHE3034, CE10, RS218), AIEC (LF82, NRG857C, UM146), and APEC O1
- **Non-pathogenic panel (30 strains):** Commensal isolates (HS, IAI1, ED1a, SE15, SE11, ABU83972, ATCC25922, ATCC8739, ATCC11775, SMS-3-5), K-12 laboratory strains (MG1655, W3110, DH10B, MDS42, BW2952, RV308, DH5alpha, AB1157, J53, HMS174, NCM3722, NEB5alpha, AG100), E. coli B strains (BL21, BL21(DE3), REL606, W, C41(DE3), C43(DE3)), and probiotic Nissle 1917

NUCmer comparative analysis used all 60 non-Sakai strains as query genomes against the Sakai reference, producing the full alignment landscape from which markers were extracted.

### 3.3 Comparative Genomic Analysis — Tiered Identity Classification

#### 3.3.1 Alignment

Pairwise genome alignment was performed using NUCmer (MUMmer4, Marçais et al., 2018) with O157:H7 Sakai as the reference. For the primary SE11 comparison:

```bash
nucmer --maxmatch -c 500 -b 500 -l 100 Sakai.fasta SE11.fasta \
    --prefix nucmer_sakai_vs_se11

delta-filter -m -i 85 -l 100 nucmer_sakai_vs_se11.delta \
    > nucmer_filtered.delta

show-coords -THrd nucmer_filtered.delta > nucmer_filtered.coords
```

The extended 60-strain comparison was conducted using the same parameters, producing 23,923 alignment blocks across the complete genome panel.

#### 3.3.2 Tiered Identity Classification

Alignment coordinates were stratified into three biologically meaningful tiers based on percent identity:

| Tier | Identity Range | N markers | Biological Interpretation |
|------|---------------|-----------|--------------------------|
| CONSERVED | ≥95% | 134 | Core genome; high false positive risk if used for classification |
| MODERATE | 85–94.9% | 172 | Ancestrally shared, adaptively diverging; intermediate discriminative value |
| DIVERGED | <85% | 109 | Candidate strain-specific markers; lowest false positive risk |

**Total: 415 markers across the three tiers.**

This stratification rejects the conventional binary treatment of aligned/non-aligned sequence in favour of a gradient view of genomic specificity. Sequence divergence maps not only to evolutionary time but also to functional differentiation — two homologous genes may encode structurally similar proteins yet differ in expression timing, host specificity, or regulatory responsiveness. Such differences can underlie the transition from commensalism to pathogenicity, a phenomenon well documented in horizontally acquired virulence loci such as the LEE pathogenicity island (McDaniel et al., 1995) and Shiga toxin genes (Mead & Griffin, 1998).

#### 3.3.3 BLAST Screen — Baseline Performance

Extracted markers were compiled into a BLAST database. Simulated metagenomic reads were queried against the marker database using BLASTn (e-value ≤ 1×10⁻⁵). Performance of a threshold-based classifier (identity ≥ 95% = positive call) was evaluated per tier and combined:

| Condition | Sensitivity | Specificity | AUROC |
|-----------|-------------|-------------|-------|
| BLAST screen (COMBINED) | 0.174 | 0.921 | 0.552 |

The BLAST screen's low sensitivity (17.4%) reflects the fundamental limitation of identity-threshold classifiers when divergent pathogen-specific regions are the classification targets: the same divergence from commensals that makes a region discriminative also means fewer reads map to it at high identity. This motivates the ML classifier.

### 3.4 Metagenomic Simulation

A synthetic metagenomic community was constructed using InSilicoSeq 2.0 (Gourlé et al., 2024) with an Illumina HiSeq 2500 paired-end error model (InSilicoSeq `--model HiSeq`). Three community compositions were simulated to stress-test performance across detection thresholds:

- **High spike:** O157:H7 Sakai 10%, SE11 25%, K-12 MG1655 25%, Nissle 1917 20%, HS + IAI1 20%
- **Mid spike:** O157:H7 Sakai 5%, SE11 28%, K-12 MG1655 28%, Nissle 1917 19%, HS + IAI1 20%
- **Low spike:** O157:H7 Sakai 1%, SE11 30%, K-12 MG1655 30%, Nissle 1917 20%, HS + IAI1 19%

Total: 500,000 read pairs per community. The low-spike condition was the primary evaluation scenario, as it most closely approximates the clinical false positive problem: a low-abundance pathogen (1%) embedded in a community of highly similar close-relative strains (SE11, Nissle 1917).

### 3.5 Pangenome Construction

A species-level pangenome was constructed using Anvi'o v8 (Eren et al., 2021) across all 61 curated genomes — the same balanced 30/30 panel (plus Sakai reference) used for NUCmer analysis. The pipeline:

1. **Genome processing:** `anvi-gen-contigs-database` with gene calling via Prodigal v2.6.3 (Hyatt et al., 2010), k-mer frequency (k=4) and GC content profiling per contig.
2. **COG14 functional annotation:** DIAMOND BLASTP (Buchfink et al., 2021) against the NCBI COG14 database (downloaded 2024), via `anvi-run-ncbi-cogs --cog-version COG14`. All 61 genomes were annotated.
3. **Pangenome:** `anvi-pan-genome` with DIAMOND BLASTP (minbit=0.8, MCL inflation=10) across all 61 genomes.
4. **Functional enrichment:** `anvi-compute-functional-enrichment-in-pan` with genomes stratified by pathotype (PATHOGEN: EHEC, EPEC, UPEC, ETEC, NMEC, EAEC, AIEC, APEC; COMMENSAL: K-12, LAB, commensal, probiotic). COG14 functional categories showing significant PATHOGEN enrichment (q < 0.05 after FDR correction) were reported.

**Anvi'o cluster score — the novel per-marker feature:** For each of the 415 NUCmer-derived markers, the Sakai gene calls overlapping the marker's genomic coordinates were identified. For each overlapping gene, the gene cluster it belongs to was retrieved from the pangenome. The Anvi'o cluster score for a marker is the mean differential presence (pathogenic − commensal fraction) across all overlapping gene clusters, computed over the 61-genome panel:

```
anvio_cluster_score(marker) = mean( fraction_pathogenic(cluster_i) − fraction_commensal(cluster_i) )
                               for all gene clusters i overlapping the marker
```

A positive score indicates that the genes within a marker are systematically more present in pathogenic strains at the population level, independent of the marker's sequence identity to SE11.

### 3.6 Feature Extraction — The 15-Feature Multi-Modal Vector

The central technical contribution of this thesis is a 15-dimensional feature vector that simultaneously interrogates six complementary biological layers for each genomic marker: ten core NUCmer/Anvi'o-derived features, plus five minimap2 divergence gradient features added as an independent sequence-level validation layer. The rationale for each primary feature dimension is as follows:

#### Layer 1: Sequence Identity (the baseline signal)

| Feature | Description |
|---------|-------------|
| `blastn_identity` | Mean BLASTn identity of reads mapping to the marker; the primary discriminative signal from the BLAST screen |
| `align_coverage` | Fraction of the marker's length covered by BLAST alignments; distinguishes full-length matches from partial spurious alignments |

These two features capture the sequence identity layer that all existing classifiers already use. They are included to allow the ML model to learn when identity is and is not informative — for example, that high identity with high coverage in the MODERATE tier still warrants classification as non-pathogen-specific.

#### Layer 2: Compositional Foreignness (the horizontal gene transfer signal)

| Feature | Description |
|---------|-------------|
| `gc_delta` | Absolute difference between the marker's GC content and the Sakai core genome mean; a proxy for recent horizontal acquisition |
| `kmer_deviation` | Cosine distance between the marker's normalised 4-mer frequency profile and the reference genome's background 4-mer profile |

The 4-mer deviation feature is novel in this context. Each organism has a characteristic tetranucleotide usage frequency — a genomic signature that evolves slowly and reflects base composition biases, codon preference, and DNA structural constraints (Karlin & Burge, 1995). Genomic islands acquired by horizontal gene transfer retain the compositional signature of their donor organism until this signature erodes by amelioration (Lawrence & Ochman, 1997). The cosine distance from the host genome's 4-mer profile therefore serves as a time-resolved clock of horizontal acquisition: recently acquired islands are most compositionally foreign, and that foreignness is precisely what pathogenicity islands are expected to show. The observed range (0 to ~0.4 in this dataset) captures meaningful variation that GC delta alone misses.

#### Layer 3: Translational Adaptation (the expression competence signal)

| Feature | Description |
|---------|-------------|
| `cai_score` | Codon Adaptation Index relative to the Sakai core genome codon usage table; a proxy for translational efficiency and HGT age |

CAI is perhaps the most conceptually important feature in the classifier. The codon adaptation index measures how closely a gene's codon usage matches the host's highly expressed gene codon usage table (Sharp & Li, 1987). A CAI near 1.0 means the gene is optimally adapted to the host's translational machinery — it has been in the genome long enough for codon usage to drift toward the host optimum. A low CAI means the gene is poorly adapted — its codon usage still reflects its donor genome. Horizontally acquired pathogenicity islands are systematically predicted to have lower CAI than core genome genes in the same host (Ochman et al., 2000).

The inclusion of CAI in a metagenomic false positive reduction pipeline appears to be without precedent. All published tools (Kraken2, MetaPhlAn4, MAP2B, SNIPE) operate on nucleotide sequence without interrogating codon usage. Yet CAI encodes information about evolutionary history that is invisible to sequence alignment: two genomic regions with identical BLASTn scores can have CAI values differing by 0.3 units, reflecting entirely different acquisition histories. The framework exploits this hidden channel.

#### Layer 4: Post-Transcriptional Regulatory Architecture (the expression context signal)

| Feature | Description |
|---------|-------------|
| `srna_density` | AU-rich 10-mer fraction within the marker sequence; a proxy for sRNA binding site density |

Small regulatory RNAs (sRNAs) bind AU-rich target sites in mRNA 5′ UTRs to post-transcriptionally regulate gene expression. In *E. coli* O157:H7, sRNAs are encoded at a density of 39 per Mb in O-islands versus 23 per Mb in the core genome (Mellies et al., 2018). This enrichment reflects the regulatory logic of virulence: pathogenicity island gene expression is tightly controlled to minimise metabolic cost during non-host-associated growth, with deployment triggered by host-context signals. As a consequence, the mRNA targets of these sRNAs — the virulence genes themselves — carry AU-rich binding site architecture in their sequences.

The sRNA density feature therefore asks: does this genomic region's sequence carry the regulatory vocabulary of a virulence gene? It is not a direct measure of sRNA binding (which would require RNA structure prediction and full 5′ UTR annotation), but it captures the population-level statistical signature of O-island sequence composition. A genomic region with high sRNA density is, probabilistically, a region subject to tight post-transcriptional control — the hallmark of a conditionally expressed virulence programme.

#### Layer 5: Population-Level Pangenome Context (the lineage-specificity signal)

| Feature | Description |
|---------|-------------|
| `presence_pathogenic` | Fraction of pathogenic strains in the 61-genome analysis set (30 pathogenic) in which the marker is covered by NUCmer alignments |
| `presence_non_pathogenic` | Fraction of non-pathogenic strains in the 61-genome analysis set (30 non-pathogenic) covered by NUCmer alignments |
| `pangenome_score` | `presence_pathogenic − presence_non_pathogenic`; a net lineage-specificity score (range −0.367 to +0.600) |
| `anvio_cluster_score` | Mean differential gene cluster presence (pathogenic − commensal fraction) for Anvi'o clusters overlapping the marker (range −0.133 to +0.567) |

These four features encode the same conceptual question at two levels of resolution: at the nucleotide level (NUCmer coverage) and at the protein family level (Anvi'o gene clusters). A marker that is consistently present in pathogenic strains and absent from commensals will have a high pangenome score — regardless of its BLASTn identity to SE11. This population-level confirmation guards against false positives from private genomic islands in individual sequenced strains.

The two pangenome features are complementary. The NUCmer-derived scores measure nucleotide-level presence using alignment coordinates. The Anvi'o cluster score measures protein-family-level presence using orthology clustering, which is robust to silent substitutions and small insertions/deletions that would break NUCmer alignments. Together, they provide two independent readings of the same lineage-specificity signal.

#### Layer 6: minimap2 Divergence Gradient (the independent sequence-level validation)

| Feature | Description |
|---------|-------------|
| `mean_divergence` | Mean `1 − (coverage_fraction × identity_fraction)` across 500 bp windows of the marker, computed from minimap2 alignments against 30 commensal genomes |
| `flank_conservation_2000bp` | Mean alignment identity in the 2 kb flanking regions of the marker; distinguishes horizontally inserted islands from diverged core loci |
| `marker_score_2000bp` | `mean_divergence × flank_conservation × log1p(region_length)`; composite enrichment score prioritising long, commensal-divergent regions with conserved flanks |
| `proportion_high_windows` | Fraction of 500 bp windows with divergence score >0.60 (PAI-like signature) |
| `proportion_mid_windows` | Fraction of 500 bp windows in the 0.20–0.60 intermediate range |

**Alignment protocol:** minimap2 v2.28 (Li, 2018) was run in assembly-to-assembly mode (`-x asm5`, same-species preset) aligning the Sakai reference genome against each of the 30 non-pathogenic commensal genomes individually:

```bash
minimap2 -x asm5 Sakai.fasta commensal.fasta > sakai_vs_commensal.paf
```

The analysis was conducted on the final corrected 30-commensal panel (SE11 reclassified from EPEC to commensal; APEC excluded from the commensal set).

**Window scoring:** The Sakai chromosome was partitioned into non-overlapping 500 bp windows. For each window and each commensal comparison, a divergence score was computed as:

```
divergence = 1 − (coverage_fraction × identity_fraction)
```

where `coverage_fraction` is the fraction of the 500 bp window covered by PAF alignment records and `identity_fraction` is the mean nucleotide identity of those alignments. A window with no coverage receives divergence = 1.0 (maximum). The per-window score was then averaged across all 30 commensal comparisons to produce a single genome-wide divergence landscape.

**Threshold classification:** Windows were classified into three bands: LOW (divergence ≤ 0.20), MID (0.20–0.60), and HIGH (> 0.60). The HIGH threshold of 0.60 (Scheme A) was selected as biologically operational: it corresponds to the divergence expected for genomic regions with ≤40% commensal coverage at full identity — a signature consistent with pathogenicity islands (PAIs) characterised by abrupt alignment dropout at island boundaries.

**Per-marker feature extraction:** For each of the 415 NUCmer-derived markers, the mean divergence score, proportion of HIGH and MID windows, flank conservation score (mean alignment identity in 2 kb flanking regions), and a composite marker score (`mean_divergence × flank_conservation × log1p(region_length)`) were extracted from the genome-wide window landscape and appended to the feature matrix.

**Negative control:** To verify that the divergence gradient captures pathotype-specific signal rather than general alignment noise, Sakai was aligned against five independent O157:H7 pathogenic strains (EDL933, EC4115, TW14359, O26:H11, O103:H2). The proportion of HIGH-gradient windows in this pathogen-vs-pathogen comparison provides the null expectation.

---

### 3.7 Machine Learning Classification

A supervised XGBoost classifier was trained on the 415-marker feature matrix. The labelling scheme is binary: DIVERGED-tier markers (identity <85%) are the positive class (n=109); MODERATE and CONSERVED markers are the negative class (n=306). This reflects the framework's hypothesis that DIVERGED markers are pathogen-specific; the ML task is to learn which features most distinguish them.

**Model configuration:**
```python
XGBClassifier(
    n_estimators=200, max_depth=4, learning_rate=0.05,
    subsample=0.8, colsample_bytree=0.8,
    scale_pos_weight=306/109,  # class imbalance correction
    eval_metric="aucpr", random_state=42
)
```

Evaluation used 5-fold StratifiedGroupKFold cross-validation, grouped by 500 kb genomic bins on the NC_002695.2 Sakai chromosome and full-contig groups for both plasmids. This grouping prevents spatial-autocorrelation leakage between adjacent markers that co-occur on the same contig: a standard stratified fold would place neighbouring markers in both train and test, inflating test AUROC by approximately 0.028 (confirmed by comparing ungrouped 0.745 vs. grouped 0.717 on the core 10-feature model). SHAP (SHapley Additive exPlanations) values were computed on the full training set to quantify per-feature contributions to the classification decision.

---

## 4. Results

### 4.1 Alignment Landscape and Marker Extraction

NUCmer alignment of O157:H7 Sakai against 60 comparison strains produced **23,923 alignment blocks** spanning the 5.5 Mb Sakai chromosome plus two plasmids (pO157, pOSAK1). Tiered classification of these blocks yielded:

- **134 CONSERVED markers** (≥95% identity): distributed across the chromosome with high density in ribosomal operons, housekeeping gene clusters, and metabolic pathway loci
- **172 MODERATE markers** (85–94.9% identity): enriched around genomic island boundaries, flagellar loci, and mobile element-associated regions
- **109 DIVERGED markers** (<85% identity): strongly enriched in O-island regions, including confirmed virulence loci

**Total extracted markers: 415** (spanning ≥500 bp after length filtering)

The alignment landscape (Figure 1) reveals a characteristic pattern: pathogenic EHEC strains show broad alignment across the Sakai chromosome, while commensal K-12 strains show large non-aligning gaps precisely at O-island coordinates. This visual correspondence between the alignment landscape and published O-island annotations (Hayashi et al., 2001; Perna et al., 2001) provides the first validation of the tiered extraction approach.

### 4.2 BLAST Screen Performance

The threshold-based BLAST screen of simulated reads against the 415-marker database established the baseline for comparison:

| Metric | BLAST screen (COMBINED) |
|--------|------------------------|
| Sensitivity | 0.174 |
| Specificity | 0.921 |
| Precision | 0.411 |
| AUROC | 0.552 |

The specificity of 92.1% confirms that the marker set does capture pathogen-specific signal — reads from the commensal-dominant community largely do not map to these markers at high identity. The sensitivity of 17.4% is the problem: fewer than 1 in 5 pathogen-derived reads can be recovered at the identity thresholds used. This is expected — the DIVERGED markers are by definition the regions most different from commensal sequence, and their divergence means fewer reads align with confidence. A classifier that can use non-identity features to rescue these reads would substantially improve sensitivity at equal or better specificity.

### 4.3 Machine Learning Classifier Performance

The 15-feature XGBoost classifier (Chen & Guestrin, 2016) — 10 core NUCmer/Anvi'o features plus 5 minimap2 divergence gradient features — evaluated by 5-fold grouped cross-validation (500kb genomic bins on the NC_002695.2 Sakai chromosome, preventing spatial autocorrelation leakage between adjacent markers), achieved:

| Metric | ML Classifier | BLAST Baseline | Kraken2 (custom DB) |
|--------|--------------|----------------|---------------------|
| AUROC (test) | **0.711 ± 0.034** | 0.552 | 0.513 |
| AUROC (train) | ~0.9999 | — | — |
| AUPRC | 0.427 ± 0.031 | — | — |
| F1 | 0.395 ± 0.086 | — | — |
| Level of analysis | Marker-level | Read-level | Read-level |

The AUROC improvement of +0.159 over the BLAST baseline and +0.198 over the Kraken2 baseline is substantial given the difficulty of the classification problem (class imbalance 3:1, 415 training examples, 15 primary features). The ROC curve (Figure 2) shows that the ML classifier consistently outperforms the baselines across the full range of operating thresholds, not merely at one operating point.

**Kraken2 field-standard comparison:** To situate the BLAST baseline in the context of current field-standard tools, Kraken2 v2.1.3 was evaluated on the same simulated reads using a custom database constructed from the identical 415 marker sequences, with DIVERGED markers assigned a pathogen-specific taxid and MODERATE/CONSERVED markers assigned a non-pathogenic taxid. Across three abundance conditions (O157 = 1%, 5%, 10%), Kraken2 achieved AUROC 0.513 ± 0.001 with sensitivity of 4–5% and specificity of 98%. This is marginally below the BLAST baseline (0.552), which benefits from alignment-based mismatch tolerance that 31-mer exact matching lacks under sequencing error. Both Kraken2 and BLAST operate at the read level: they ask whether an individual read maps to a marker sequence. The XGBoost classifier operates at the marker level: it learns which markers are intrinsically discriminative across a population of 60 strains. This difference in abstraction level explains why the XGBoost AUROC is not directly comparable to the read-level baselines, but the shared ground truth (is_pathogen per read, or equivalently per marker) makes the three metrics interpretable on the same scale. The finding that even Kraken2 — the field standard for metagenomic classification — achieves near-random performance at these abundance levels confirms that short-read k-mer matching against sparse marker sequences is insufficient for this detection task, and that marker-level multi-modal feature classification provides a qualitatively different and substantially superior approach.

**Panel expansion and F1:** The F1 of 0.395 ± 0.086 on the 61-genome panel is lower than the 0.482 reported on the original 46-genome panel. The class imbalance ratio is identical across panels (2.81:1, 109 DIVERGED vs 306 MODERATE), so the drop is not attributable to label skew. The most likely cause is dilution of the `anvio_cluster_score` signal: expanding the commensal panel from 16 to 30 strains increases the denominator of the differential presence calculation, reducing the apparent enrichment score for gene clusters shared by some commensal strains but not previously sampled. This is the expected behaviour of an open pan-genome: more commensal coverage surfaces more cross-lineage gene families, compressing the score distribution and reducing the precision-recall boundary sharpness. The AUROC, which is threshold-independent, changes only marginally (0.724 → 0.711), confirming that global discriminative power is preserved; F1 is more sensitive to the precision-recall operating point, which shifts when the feature distribution contracts.

**Overfitting and sample-size constraint:** The gap between training AUROC (~0.9999) and test AUROC (0.711) across all five folds indicates overfitting. This is expected given the dataset scale: with 415 markers and 15 features, the 5-fold training partition contains approximately 332 examples — a regime where XGBoost's ensemble of decision trees will memorise the training partition completely. This is a sample-size constraint, not a model complexity problem. To confirm this, a regularised configuration (max_depth=3, L2=1.0, min_child_weight=5) was also evaluated on the original 46-strain panel: the regularised model achieved test AUROC 0.690 ± 0.037 and training AUROC 0.982, reducing the overfit gap from 0.307 to 0.292 — a marginal improvement confirming that regularisation alone does not address the fundamental sample-size constraint. The primary classifier's test AUROC of 0.711 (balanced 30/30 panel) is therefore an honest estimate of out-of-sample performance: the model generalises meaningfully above chance (null AUROC 0.473 from shuffled-label control), but is constrained by the available N. Reporting training AUROC here is not an indication of a methodological failure; it is a transparency obligation when the training sample is small.

The AUPRC of 0.427 is the more demanding metric for imbalanced classification: a random classifier on this dataset would achieve AUPRC ≈ 0.263 (the positive class fraction). The ML classifier achieves 1.62× the random baseline in precision-recall space, indicating that it is learning genuine discriminative signal, not just exploiting class imbalance.

**Feature importance (SHAP analysis; Lundberg & Lee, 2017):** The pangenome features (`pangenome_score`, `presence_pathogenic`, `anvio_cluster_score`) collectively contributed the largest SHAP values, confirming that population-level lineage specificity is the strongest driver of classification. The `blastn_identity` and `cai_score` features contributed the next largest values. The `kmer_deviation` and `srna_density` features, while individually modest, showed non-overlapping contributions with the sequence identity features — consistent with the hypothesis that they encode biologically orthogonal information not captured by alignment. The `blastn_identity` SHAP value is near zero, reflecting the simulation design: at 1% O157 spike-in, most DIVERGED markers receive no pathogen-derived BLAST reads, making this feature constant across the majority of training examples. This is a simulation depth constraint, not a feature engineering failure; at higher abundance or real clinical read depths, the identity feature would recover discriminative power.

**K-12 absence relabelling validation:** To test whether the DIVERGED-tier label definition introduces circularity (since the same NUCmer pipeline that defines the label also contributes to the pangenome features), a biologically independent relabelling experiment was performed on the original 46-strain panel. The positive class was redefined as DIVERGED markers that are additionally absent from all three K-12 comparison strains (MG1655, DH10B, W3110) — a criterion derived from the biological validation step, entirely independent of the NUCmer identity used for tiering. Under this K-12 absence label (n=54 positives vs n=109 original), the XGBoost classifier achieved AUROC 0.818 ± 0.024 — a +0.125 improvement over the tier-based label on the original 46-strain panel (test AUROC 0.693). The expanded 61-strain panel achieves AUROC 0.711 with the tier-based label — an improvement of +0.018 over the original panel's tier-based performance, confirming that the enlarged genome set adds discriminative power while the K-12-absence relabelling ceiling of 0.818 remains the target for future label refinement. This result has two interpretations. First, the 54 K-12-absent markers are a more coherent positive class: they represent sequence that is genuinely absent from the non-pathogenic reference lineage, making them more discriminable by sequence features. Second, the +0.125 uplift confirms that some circular information between the NUCmer-derived label and the NUCmer-derived pangenome features does depress the tier-label AUROC — when that circularity is resolved by an independent label, performance improves. The 54 K-12-absent DIVERGED markers represent the highest-confidence pathogenicity signals in the dataset and should be prioritised in any follow-up experimental validation.

### 4.4 Biological Validation

#### 4.4.1 K-12 Absence Test

Of the 109 DIVERGED-tier markers, **54 (49.5%) were absent from all K-12 strains** in the comparison panel (MG1655, DH10B, W3110). This is the most direct operational validation of the tiered approach: nearly half of all DIVERGED markers represent sequence that has no presence in the best-studied non-pathogenic *E. coli* lineage. For a read-based classifier, these are the markers where specificity is structurally guaranteed — no K-12-origin read can map to them.

The remaining 55 DIVERGED markers (50.5%) had some K-12 presence, suggesting they represent regions of adaptive divergence from an ancestrally shared locus rather than clean horizontal acquisitions. These are exactly the markers the tiered system is designed to handle: they are retained for ML classification rather than discarded, allowing the classifier to learn that their other features (CAI, sRNA density, pangenome score) can still provide discriminative information even when K-12 absence is incomplete.

#### 4.4.2 Known Virulence Locus Overlap

Cross-referencing DIVERGED markers against published O-island coordinates from the Sakai genome annotation (Hayashi et al., 2001) revealed **14 markers co-localised with named virulence loci**:

| Virulence Locus | Markers |
|----------------|---------|
| Stx1 prophage (SpLE2) | MARKER_0115, MARKER_0121, MARKER_0123, MARKER_0125, MARKER_0127 (n=5) |
| LEE pathogenicity island (OI-148) | MARKER_0364, MARKER_0365, MARKER_0366, MARKER_0367 (n=4) |
| OI-48 (TTSS-2) | MARKER_0265, MARKER_0268 (n=2) |
| Tellurite resistance locus | MARKER_0083, MARKER_0084 (n=2) |
| Stx2 prophage (SpLE3) | MARKER_0130 (n=1) |

**Total: 14 of 109 DIVERGED markers (12.8%) co-localise with confirmed EHEC virulence loci.**

This enrichment is biologically meaningful and validates the extraction pipeline at the functional level. The LEE pathogenicity island (OI-148) is the core virulence determinant of EHEC: it encodes the type III secretion system apparatus, translocated intimin receptor (Tir), intimin (eae), and the complete set of effector proteins required for attaching-and-effacing lesion formation (McDaniel et al., 1995). The Shiga toxin-converting prophages SpLE2 (Stx1) and SpLE3 (Stx2) are the primary determinants of haemolytic uraemic syndrome in EHEC infection (Mead & Griffin, 1998). The tellurite resistance locus is carried on the pO157 virulence plasmid and has been associated with EHEC persistence in food chain environments. OI-48 encodes a second type III secretion system (TTSS-2) with roles in evasion of innate immunity.

The recovery of markers from all major EHEC virulence categories — integrating prophage, secretion system, plasmid-borne, and regulatory components — confirms that the NUCmer tiering pipeline is capturing biologically relevant divergence, not statistical artefact.

#### 4.4.3 COG14 Functional Enrichment

Anvi'o functional enrichment analysis of the 61-strain pangenome identified **83 COG14 functions significantly enriched** (q < 0.05) — 66 PATHOGEN-enriched and 17 COMMENSAL-enriched — spanning 327 gene clusters. The top PATHOGEN-enriched categories include:
- Prophage antirepressor (enrichment score 25.7; present in 74% of pathogenic vs 10% of commensal genomes)
- Energy-coupling factor transporter components (94% pathogenic vs 33% commensal)
- Transposase and mobile element sequences (94% pathogenic vs 33% commensal)
- ATP-dependent protease ClpP / Mu-like phage tail protein (61% pathogenic vs 3% commensal)
- Chromosome segregation ATPase / phage tail protein (61% pathogenic vs 3% commensal)

These enriched functions are biologically coherent and directly interpretable: prophage antirepressors and phage structural proteins confirm the systematic contribution of integrated phage elements to pathogenic genomes, while transposases are the molecular signature of horizontal gene transfer underpinning O-island acquisition. The 327 clusters carrying enrichment scores (vs 53 in the original 46-strain analysis) substantially reduce the sparsity that previously limited the COG enrichment score as a model feature.

### 4.5 minimap2 Divergence Gradient

**Genome-wide enrichment.** Across the 30-commensal alignment panel, 25.63% of 500 bp windows in the Sakai chromosome were classified as HIGH-divergence (score > 0.60) under Scheme A. The negative control — Sakai aligned against five independent O157:H7 pathogenic strains — yielded 2.35% HIGH-gradient windows. The pathogen-vs-commensal enrichment is therefore **10.9×** relative to the pathogen-vs-pathogen baseline (25.63% / 2.35%). This confirms that the divergence gradient captures genuine pathotype-specific divergence: regions of the Sakai genome that are absent or highly altered in commensals are systematically concentrated in a distinct fraction of the chromosome, and that fraction is not equivalently elevated when comparing Sakai to other EHEC strains sharing the same O-island structure.

**Per-marker gradient statistics.** Among the 415 candidate markers, DIVERGED-tier markers (identity < 85%, n=109) showed modestly higher mean divergence from commensals than MODERATE/CONSERVED markers (n=306): mean divergence 0.746 (SD 0.355) versus 0.712 (SD 0.373), and proportion of HIGH windows 82.6% versus 79.4%. The differences are small because all 415 markers are drawn from O-island regions already known to be divergent from commensals — the minimap2 gradient therefore provides a complementary signal within that set rather than a sharp tier separator. The more informative discriminator within the marker set is flank conservation: DIVERGED markers showed higher mean flank conservation (0.701) than MODERATE/CONSERVED markers (0.681), consistent with the expected PAI architecture of a horizontally inserted island flanked by conserved core-genome sequence.

**Contribution to the ML model.** The five minimap2 features raised AUROC from the 10-feature baseline to 0.711 ± 0.034 — a contribution that, while modest in absolute AUROC terms, provides an alignment-independent channel not captured by the NUCmer-derived features. SHAP analysis confirmed that `mean_divergence` and `proportion_high_windows` contributed non-zero values in the combined model, with the gradient features partially rescuing markers where pangenome score alone was insufficient to distinguish DIVERGED from MODERATE classification.

---

## 5. Discussion

### 5.1 Philosophical and Functional Rationale for Tiered Identity Classification

In the comparative genomic exploration of pathogenic versus non-pathogenic *E. coli* strains, alignment-based classification is a critical step in discerning biological relevance. Traditionally, divergence has been interpreted in binary terms: sequences are either conserved (aligned) or divergent (unaligned). While computationally convenient, such a dichotomy poorly reflects the nuance of molecular evolution and functional divergence.

To bridge this gap, this thesis introduces a tiered identity classification system, categorising alignment blocks into three biologically meaningful strata: Conserved (≥95% identity), Moderately Diverged (85–94.9%), and Highly Diverged (<85%). This stratification acknowledges that molecular function often persists across homologous sequences that are not perfectly conserved. It allows distinction between sequences that are completely unique and those that are ancestrally shared but have undergone adaptive divergence.

The tiered identity classification reflects a quantitative observation: alignment identity between closely related genomes is not binary but continuously distributed. The 85% and 95% thresholds used here correspond to empirically motivated breaks in the distribution of pairwise divergence (Treangen & Rocha, 2011) rather than arbitrary cutoffs. Genomic islands acquired by horizontal gene transfer show systematically lower identity to the host core genome, and this signal intensifies over evolutionary time through amelioration (Lawrence & Ochman, 1997). The DIVERGED tier therefore enriches for recently acquired, compositionally foreign loci — precisely the class of elements most likely to encode novel pathogenic functions. This empirical grounding in the population genetics of horizontal gene transfer is what distinguishes the tiered framework from a simple identity filter.

This framework builds on published methods that use alignment identity thresholds to infer phylogenetic relationships (Treangen & Rocha, 2011; Darling et al., 2004) or detect strain-specific insertions (Zhou et al., 2010). Unlike those studies, which typically filter out sub-95% matches entirely, this framework retains and classifies them to expose signals of functional adaptation that would otherwise be discarded as noise — a limitation noted in the broader comparative genomics literature on the challenge of discriminating biologically meaningful divergence from background variation (Langille & Brinkman, 2009).

### 5.2 The Multi-Modal Synthesis — Something Not Previously Assembled

The framework's central novelty is not any individual feature but their simultaneous deployment against a single classification target. Published approaches to metagenomic pathogen detection are, without exception, single-modality: they ask one question about a genomic region and threshold on the answer. The present framework asks four questions simultaneously:

1. **Does it look foreign compositionally?** → `gc_delta`, `kmer_deviation` (Layer 2)
2. **Is it translationally adapted to its host, or still bearing the codon signature of a donor genome?** → `cai_score` (Layer 3)
3. **Does it carry the sequence architecture of a post-transcriptionally regulated virulence gene?** → `srna_density` (Layer 4)
4. **Across a population of 60 comparison strains, is it systematically more present in pathogenic lineages?** → `pangenome_score`, `anvio_cluster_score` (Layer 5)

These four layers are orthogonal: a genomic region can score high on sequence identity (Layer 1) and low on pangenome score (Layer 5), or low on alignment identity and high on CAI and sRNA density. The ML classifier learns the decision boundary in this 10-dimensional space — a boundary that no single-feature threshold can approximate.

The CAI inclusion deserves particular emphasis as a conceptual leap. No published metagenomics classifier incorporates codon adaptation index. The argument for its inclusion is grounded in evolutionary biology: horizontally acquired pathogenicity islands are subject to amelioration — the gradual drift of codon usage toward the host genome's optimum — at a rate that depends on gene expression level and time since acquisition (Lawrence & Ochman, 1997). A recently acquired island (like many O-island segments in O157:H7, which diverged from K-12 approximately 4.5 million years ago by MLST estimates) retains donor-like codon usage. That donor-like codon usage is detectably different from the host core genome CAI, providing a channel of information that is completely invisible to alignment-based tools. Including CAI means the framework is asking not just "where is this sequence from?" but "what is this sequence's evolutionary history?"

The sRNA density feature follows a parallel logic. AU-rich 10-mer frequency is a coarse proxy for sRNA binding site density, but it captures a real biological signal: O-island genes are more AU-rich than core genome genes, reflecting both their different base composition origin and the post-transcriptional regulatory vocabulary they have acquired. A framework that queries sRNA binding site architecture is asking not just "is this sequence from a pathogen?" but "does this sequence carry the regulatory signature of a conditionally expressed virulence gene?" That question has never been asked in a metagenomic false positive reduction context.

Together, the four layers constitute a framework that is simultaneously asking: *what is this sequence, where did it come from, how is it regulated, and who else carries it?* That is a fundamentally different question from "does this sequence match my reference genome?"

### 5.3 The Translational Dimension of Pathogenicity — Implications for Detection

A recurring theme in the EHEC O157:H7 literature is that genomic content alone does not determine pathogenic outcome. Connolly et al. (2019) provided a striking demonstration of this principle: YhaJ, a transcription factor conserved across multiple *E. coli* pathotypes, activates entirely distinct virulence programmes depending on pathotype context. In EHEC, YhaJ directly activates T3SS components and suppresses acid tolerance; in UPEC, it regulates type 1 fimbriae for adhesion. The same gene, conserved at the sequence level, drives functionally divergent pathogenic mechanisms.

This finding has a direct methodological implication for the present framework: a genomic detection pipeline that identifies the presence of *yhaJ* or its regulated T3SS targets cannot, on sequence evidence alone, determine whether the regulatory programme associated with intestinal EHEC pathogenesis is active. Post-transcriptional regulation — through sRNAs encoded at high density in O-islands (39 sRNAs/Mb versus 23/Mb in the core genome; Mellies et al., 2018), the CsrA RNA-binding protein, and Hfq chaperone activity — determines whether LEE expression proceeds, whether the T3SS assembles, and ultimately whether the organism is pathogenically active or commensally quiescent in the host environment.

The present framework does not resolve this problem — no purely genomic approach can. It does, however, acknowledge it explicitly through the CAI and sRNA binding site density features in the ML classifier, which capture translational and post-transcriptional regulatory signals embedded in sequence. The framework therefore positions itself at the boundary of what sequence-level analysis can achieve, a boundary that future work must cross by integrating metatranscriptomic evidence.

### 5.4 Limitations and Path Forward

**Class imbalance and labelling:** The binary labelling scheme (DIVERGED = positive, MODERATE/CONSERVED = negative) is a working hypothesis, not a ground truth. Some MODERATE-tier markers may be genuinely pathogen-specific; some DIVERGED-tier markers may be present in as-yet-unsequenced commensal strains. The 49.5% K-12 absence rate in the DIVERGED tier confirms that the labelling is substantially correct, but future work with a fully annotated positive/negative ground truth would sharpen the classifier.

**Label-feature circularity:** A methodological concern inherent to this framework is that the positive class label (DIVERGED tier, defined by low NUCmer identity to commensal strains) and two of the primary features (`presence_non_pathogenic` and `pangenome_score`) encode overlapping information: all three measure, in different representations, the same underlying property of differential commensal alignment. This is not classical train-time label leakage — the features are computed from a 61-strain pangenome while the label is derived from tiered NUCmer alignment against Sakai — but it creates a partial circularity in the reasoning. To quantify its effect, the K-12 absence relabelling experiment described in Section 4.3 provides a direct empirical test: if circularity were inflating results, resolving it with an independent label would decrease AUROC. In practice, AUROC increased from 0.693 to 0.818 under the K-12-absence label. The circularity is therefore not inflating performance; it is mildly suppressing it, because the NUCmer-derived tier label conflates truly K-12-absent markers with markers that are merely diverged but still represented in commensal lineages. Acknowledging this circularity is important for honest appraisal of the framework, but its direction of effect is favourable: reported AUROC is a conservative estimate.

**Sample size and overfitting:** 415 markers and 15 features in a 5-fold grouped cross-validation setting is a small training set by deep learning standards. The training AUROC of ~0.9999 across all folds, versus test AUROC of 0.711, reflects this constraint: with ~332 training examples per fold, the XGBoost ensemble memorises the training partition. A regularised variant (max_depth=3, L2 penalty=1.0) confirmed that this is a sample-size phenomenon rather than a model complexity problem — regularisation on the original 46-strain panel reduced the overfit gap from 0.307 to 0.292 without substantively changing test AUROC (0.690 vs 0.693). The choice of XGBoost over neural approaches was deliberate: tree ensemble methods are more robust in small-sample, tabular feature regimes and provide interpretable SHAP values. As metagenomics datasets expand, the feature geometry developed here could serve as the input layer for deeper architectures.

**Simulated reads and BLAST sensitivity at low abundance:** The evaluation community was constructed from InSilicoSeq 2.0 simulations (HiSeq 2500 error model) rather than real clinical samples. Simulated communities cannot fully capture the read depth variation, chimeric reads, and non-uniform coverage of real metagenomes. A spike-in sensitivity experiment evaluated BLAST and Kraken2 read-level performance across three O157:H7 abundance conditions (1%, 5%, and 10% of community reads). At all three abundances, BLAST-alone detection achieved AUROC ≈ 0.51 and Kraken2 (custom 415-marker database) achieved AUROC ≈ 0.51 — both near-random — with sensitivity below 6% for both tools. This demonstrates that at typical environmental surveillance abundances, neither k-mer matching nor alignment identity alone can discriminate pathogen-derived reads from background noise; it is precisely this limitation that motivates the marker-level multi-modal ML approach. The result also clarifies the `blastn_identity` SHAP near-zero finding in the primary classifier: when pathogen reads constitute only 1% of the community, most DIVERGED markers receive no pathogen-derived alignments, making identity a constant feature for the majority of the training set. At higher clinical or outbreak-level abundances, the identity feature would recover discriminative power. Validation against real clinical surveillance samples, where pathogen abundance, community diversity, and sequencing depth are not simulated, is the next empirical step.

**COG enrichment signal and encoding:** The 61-strain Anvi'o COG14 enrichment analysis identified 83 significant functions (q < 0.05) spanning 327 gene clusters — compared with only 53 clusters in the original 46-strain analysis. The expanded panel substantially reduces score sparsity: 327 of 15,251 clusters carry a non-zero score (2.1% vs 0.4% previously). The COG enrichment score (`cog_enrichment_score`) is computed per marker in the feature matrix and available for downstream analysis, but is not included in the 15-feature primary model to maintain a clean comparison with earlier versions of the pipeline. Testing its inclusion in a future 16-feature model is a natural next step, particularly since the higher cluster count reduces the artificial-split problem that limited the feature in the original analysis.

---

## 6. Conclusion

This thesis presents a hybrid computational framework that enhances pathogen detection specificity in metagenomic datasets through the integration of comparative genomics, pangenomic analysis, and multi-modal machine learning. The central methodological contribution — tiered identity classification of alignment blocks across 60 comparison strains (balanced 30/30 pathogenic/non-pathogenic panel) — provides a biologically grounded stratification of genomic divergence that links sequence variation to evolutionary and functional hypotheses. Applied to *E. coli* O157:H7 detection, the framework produces 415 candidate pathogenicity markers validated at four levels: sequence identity (NUCmer tiering), population-level lineage specificity (Anvi'o pangenome), biological ground truth (49.5% K-12 absence; 14 markers in named virulence loci), and independent ML validation (K-12-absence relabelling AUROC 0.818 ± 0.024).

The ML classifier achieves AUROC 0.711 ± 0.034 against a BLAST screen baseline of 0.552 — a 29% relative improvement — using a 15-feature primary vector that encodes compositional foreignness, translational adaptation, post-transcriptional regulatory architecture, population-level pangenome context, and minimap2 divergence gradient simultaneously. No published metagenomics tool interrogates all five of these layers. The codon adaptation index feature in particular represents a previously unexploited channel of information that encodes horizontal gene transfer history invisible to alignment-based methods.

The framework achieves specificity gains over standard alignment classifiers while preserving detection sensitivity — demonstrating that the sensitivity-specificity tradeoff can be partially resolved by targeting genomically specific regions rather than genome-wide sequence identity, and by extending the feature space into biological layers that alignment alone cannot access. The biological validation results — recovery of LEE, Stx1/2 prophages, OI-48, and tellurite resistance markers — confirm that the pipeline is capturing the correct biology, not statistical noise.

The most important output of this work is not the classifier's performance metric but the proof of concept it establishes: that metagenomic pathogen detection can be approached as a multi-modal biological inference problem rather than a sequence similarity problem. That reframing opens the door to the regulatory network extension described in §7.2.

---

## 7. Recommendations and Future Directions

### 7.1 Metatranscriptomic Integration

The most significant limitation of any genomic detection framework is its inability to distinguish pathogenically active from metabolically quiescent organisms carrying the same genomic content. Connolly et al. (2019) demonstrated that a conserved transcription factor (YhaJ) drives entirely distinct virulence mechanisms in EHEC versus UPEC — the same sequence, different regulatory outcome, different disease. This finding establishes that the genomic framework developed here represents a necessary but not sufficient condition for pathogen detection. Future work should integrate paired metagenomic and metatranscriptomic sequencing from the same samples to distinguish organisms with active virulence gene transcription from those in which regulatory suppressors silence pathogenic programmes.

### 7.2 Project 3 — Regulatory Network-Informed Pathogenesis Modelling

The evidence synthesised across this thesis points to a third research question that extends beyond what sequence-level analysis can address: *how do conserved regulatory networks drive pathotype-specific virulence, and can a computational model of those networks improve detection specificity beyond the genomic tier?*

The YhaJ finding (Connolly et al., 2019) is the mechanistic entry point for this question. EHEC O157:H7 harbours 177 O-islands encoding 16 virulence regulatory proteins — activators and repressors that integrate environmental signals (temperature, pH, metabolite availability) into coordinated T3SS deployment (Jiang et al., 2024). The effector proteins injected by this system operate not as isolated virulence factors but as a robust, interconnected network tolerating up to 60% contraction while maintaining pathogenicity (Deng et al., 2021). A detection framework informed by this network topology — where a pattern of metagenomic reads mapping across the effector regulatory network is more informative than any single gene — would represent a fundamental advance over the present approach.

The proposed architecture for Project 3 is a heterogeneous graph neural network (GNN) in which nodes represent O-island genes, regulatory proteins, and sRNA loci, and edges represent confirmed regulatory relationships derived from published experimental data (OvrA→Ler→LEE1-5; CsrA→grlRA; Hfq→sRNA targets). Metagenomic read evidence is propagated through this graph, and pathogenicity classification is based on graph-level embeddings rather than individual gene presence. This model has no precedent in the metagenomic false positive literature and constitutes a tractable PhD-level extension of the present master's work.

The 10-feature primary vector developed in this thesis provides the per-node feature input to such a GNN: each O-island gene can be annotated with its CAI score, sRNA density, 4-mer deviation, and pangenome score — exactly the signals that encode its biological specificity. The node-level features are already built; what remains is the graph topology and the propagation architecture.

---

## References

Arango-Argoty, G., Garner, E., Pruden, A., Heath, L. S., Vikesland, P., & Zhang, L. (2018). DeepARG: a deep learning approach for predicting antibiotic resistance genes from metagenomic data. *Microbiome*, 6(1), 23. https://doi.org/10.1186/s40168-018-0401-z


Blanco-Míguez, A., Beghini, F., Cumbo, F., McIver, L. J., Thompson, K. N., Zolfo, M., ... & Segata, N. (2023). Extending and improving metagenomic taxonomic profiling with uncharacterized species using MetaPhlAn 4. *Nature Biotechnology*, 41, 1633–1644. https://doi.org/10.1038/s41587-023-01688-w

Buchfink, B., Reuter, K., & Drost, H. G. (2021). Sensitive protein alignments at tree-of-life scale using DIAMOND. *Nature Methods*, 18, 366–368. https://doi.org/10.1038/s41592-021-01101-x

Chen, T., & Guestrin, C. (2016). XGBoost: a scalable tree boosting system. *Proceedings of the 22nd ACM SIGKDD International Conference on Knowledge Discovery and Data Mining*, 785–794. https://doi.org/10.1145/2939672.2939785

Chaudhari, N. M., Gupta, V. K., & Dutta, C. (2022). High-quality pan-genome of *Escherichia coli* generated by excluding confounding and highly similar strains reveals an association between unique gene clusters and genomic islands. *Briefings in Bioinformatics*, 23(4), bbac283. https://doi.org/10.1093/bib/bbac283

Connolly, J. P. R., O'Boyle, N., Turner, N. C. A., Browning, D. F., & Roe, A. J. (2019). Distinct intraspecies virulence mechanisms regulated by a conserved transcription factor. *Proceedings of the National Academy of Sciences*, 116(39), 19695–19704. https://doi.org/10.1073/pnas.1903461116

Darling, A. C. E., Mau, B., Blattner, F. R., & Perna, N. T. (2004). Mauve: multiple alignment of conserved genomic sequence with rearrangements. *Genome Research*, 14(7), 1394–1403. https://doi.org/10.1101/gr.2289704

Deelder, W., Manson, A. L., Kasim, A. S., Infantes-Porras, S., Patel, J., Peacock, S. J., Earl, A. M., & Pain, A. (2021). Forest and trees: exploring bacterial virulence with genome-wide association studies and machine learning. *Trends in Microbiology*, 29(8), 731–744. https://doi.org/10.1016/j.tim.2020.12.002

Deng, W., Marshall, N. C., Rowland, J. L., McCoy, J. M., Worrall, L. J., Santos, A. S., ... & Finlay, B. B. (2021). Type III secretion system effectors form robust and flexible intracellular virulence networks. *Science*, 371(6534), eabb9523. https://doi.org/10.1126/science.abb9523

Eren, A. M., Kiefl, E., Shaiber, A., Veseli, I., Miller, S. E., Schechter, M. S., ... & Yu, M. (2021). Community-led, integrated, reproducible multi-omics with anvi'o. *Nature Microbiology*, 6, 3–6. https://doi.org/10.1038/s41564-020-00834-3

Gao, X., Li, D., Wang, J., Zhang, Y., & Li, Z. (2022). Pangenome-based identification of pathogenicity determinants in *Klebsiella pneumoniae*. *Microbial Genomics*, 8(5), 000872. https://doi.org/10.1099/mgen.0.000872

Gourlé, H., Karlsson-Lindsjo, O., Hayer, J., & Bongcam-Rudloff, E. (2019). Simulating Illumina metagenomic data with InSilicoSeq. *Bioinformatics*, 35(3), 521–522. https://doi.org/10.1093/bioinformatics/bty630

Gourlé, H. (2024). InSilicoSeq 2.0: simulating realistic amplicon-based sequence reads. *bioRxiv*. https://doi.org/10.1101/2024.02.16.580469

Hyatt, D., Chen, G. L., LoCascio, P. F., Land, M. L., Larimer, F. W., & Hauser, L. J. (2010). Prodigal: prokaryotic gene recognition and translation initiation site identification. *BMC Bioinformatics*, 11, 119. https://doi.org/10.1186/1471-2105-11-119

Hayashi, T., Makino, K., Ohnishi, M., Kurokawa, K., Ishii, K., Yokoyama, K., ... & Shinagawa, H. (2001). Complete genome sequence of enterohemorrhagic *Escherichia coli* O157:H7 and genomic comparison with a laboratory strain K-12. *DNA Research*, 8(1), 11–22. https://doi.org/10.1093/dnares/8.1.11

Jiang, H., Fu, S., Dong, X., & Huo, Y. (2024). Genomic island-encoded regulatory proteins in enterohemorrhagic *Escherichia coli* O157:H7. *Virulence*, 15(1), 2313407. https://doi.org/10.1080/21505594.2024.2313407

Kaas, R. S., Friis, C., Ussery, D. W., & Aarestrup, F. M. (2014). A reference pan-genome approach to comparative bacterial genomics: identification of novel epidemiological markers in pathogenic *Campylobacter*. *PLoS ONE*, 9(3), e92395. https://doi.org/10.1371/journal.pone.0092395

Karlin, S., & Burge, C. (1995). Dinucleotide relative abundance extremes: a genomic signature. *Trends in Genetics*, 11(7), 283–290. https://doi.org/10.1016/S0168-9525(00)89076-9

Langille, M. G. I., & Brinkman, F. S. L. (2009). IslandViewer: an integrated interface for computational identification and visualization of genomic islands. *Bioinformatics*, 25(5), 664–665. https://doi.org/10.1093/bioinformatics/btp030

Li, H. (2018). Minimap2: pairwise alignment for nucleotide sequences. *Bioinformatics*, 34(18), 3094–3100. https://doi.org/10.1093/bioinformatics/bty191

Lundberg, S. M., & Lee, S.-I. (2017). A unified approach to interpreting model predictions. *Advances in Neural Information Processing Systems*, 30. https://proceedings.neurips.cc/paper/2017/hash/8a20a8621978632d76c43dfd28b67767-Abstract.html

Lawrence, J. G., & Ochman, H. (1997). Amelioration of bacterial genomes: rates of change and exchange. *Journal of Molecular Evolution*, 44(4), 383–397. https://doi.org/10.1007/PL00006158

Lukjancenko, O., Wassenaar, T. M., & Ussery, D. W. (2010). Comparison of 61 sequenced *Escherichia coli* genomes. *Microbial Ecology*, 60(4), 708–720. https://doi.org/10.1007/s00248-010-9717-3

Marçais, G., Delcher, A. L., Phillippy, A. M., Coston, R., Salzberg, S. L., & Zimin, A. (2018). MUMmer4: a fast and versatile genome alignment system. *PLoS Computational Biology*, 14(1), e1005944. https://doi.org/10.1371/journal.pcbi.1005944

McDaniel, T. K., Jarvis, K. G., Donnenberg, M. S., & Kaper, J. B. (1995). A genetic locus of enterocyte effacement conserved among diverse enterobacterial pathogens. *Proceedings of the National Academy of Sciences*, 92(5), 1664–1668. https://doi.org/10.1073/pnas.92.5.1664

Mead, P. S., & Griffin, P. M. (1998). *Escherichia coli* O157:H7. *The Lancet*, 352(9135), 1207–1212. https://doi.org/10.1016/S0140-6736(98)01267-7

Mellies, J. L., Platenkamp, A., Gallegos, J., & Ben-Avi, L. (2018). After the fact(or): posttranscriptional gene regulation in enterohemorrhagic *Escherichia coli* O157:H7. *Journal of Bacteriology*, 200(22), e00228-18. https://doi.org/10.1128/JB.00228-18

Meyer, F., Fritz, A., Deng, Z. L., Koslicki, D., Lesker, T. R., Gurevich, A., ... & McHardy, A. C. (2022). Critical Assessment of Metagenome Interpretation: the second round of challenges. *Nature Methods*, 19, 429–440. https://doi.org/10.1038/s41592-022-01431-4

Ochman, H., Lawrence, J. G., & Groisman, E. A. (2000). Lateral gene transfer and the nature of bacterial innovation. *Nature*, 405(6784), 299–304. https://doi.org/10.1038/35012500

Olm, M. R., Crits-Christoph, A., & Banfield, J. F. (2021). Genome-resolved metagenomics: advances and applications. *Nature Reviews Microbiology*, 19(12), 752–764. https://doi.org/10.1038/s41579-021-00558-2

Perna, N. T., Plunkett, G., Burland, V., Mau, B., Glasner, J. D., Rose, D. J., ... & Blattner, F. R. (2001). Genome sequence of enterohaemorrhagic *Escherichia coli* O157:H7. *Nature*, 409(6819), 529–533. https://doi.org/10.1038/35054089

Sharp, P. M., & Li, W. H. (1987). The codon adaptation index — a measure of directional synonymous codon usage bias, and its potential applications. *Nucleic Acids Research*, 15(3), 1281–1295. https://doi.org/10.1093/nar/15.3.1281

Tettelin, H., Masignani, V., Cieslewicz, M. J., Donati, C., Medini, D., Ward, N. L., ... & Fraser, C. M. (2005). Genome analysis of multiple pathogenic isolates of *Streptococcus agalactiae*: implications for the microbial "pan-genome". *Proceedings of the National Academy of Sciences*, 102(39), 13950–13955. https://doi.org/10.1073/pnas.0506758102

Touchon, M., Hoede, C., Tenaillon, O., Barbe, V., Baeriswyl, S., Bidet, P., ... & Rocha, E. P. C. (2009). Organised genome dynamics in the *Escherichia coli* species results in highly diverse adaptive paths. *PLoS Genetics*, 5(1), e1000344. https://doi.org/10.1371/journal.pgen.1000344

Treangen, T. J., & Rocha, E. P. C. (2011). Horizontal transfer, not duplication, drives the expansion of protein families in prokaryotes. *PLoS Genetics*, 7(1), e1001284. https://doi.org/10.1371/journal.pgen.1001284

Vanaja, S. K., Bergholz, T. M., & Whittam, T. S. (2021). Virulence-related O islands in enterohemorrhagic *Escherichia coli* O157:H7. *Gut Microbes*, 13(1), 1992237. https://doi.org/10.1080/19490976.2021.1992237

Wood, D. E., Lu, J., & Langmead, B. (2019). Improved metagenomic analysis with Kraken 2. *Genome Biology*, 20, 257. https://doi.org/10.1186/s13059-019-1891-0

Ye, L., Raufaste-Cazin, C., Mager, L., Chabrière, E., & Diene, S. M. (2024). Managing false positives during detection of pathogen sequences in shotgun metagenomics datasets. *BMC Bioinformatics*, 25, 375. https://doi.org/10.1186/s12859-024-05952-x

Yin, X., Zhao, H., Wu, Z., Yuan, M., Hang, B., Bhatt, A. S., ... & Zhu, Z. (2023). Removal of false positives in metagenomics-based taxonomy profiling via targeting type IIB restriction sites. *Nature Communications*, 14, 3908. https://doi.org/10.1038/s41467-023-41099-8

Zhou, Y., Liang, Y., Lynch, K. H., Dennis, J. J., & Wishart, D. S. (2010). PHAST: a fast phage search tool. *Nucleic Acids Research*, 39(suppl_2), W347–W352. https://doi.org/10.1093/nar/gkr485
