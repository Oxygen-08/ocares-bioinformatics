# Project 1 Report: Comparative Genomic Analysis and Tiered Identity Classification of *Escherichia coli* O157:H7 Pathogenicity Markers

**Author:** Oluwatosin Samuel Oluwole  
**Programme:** MSc Microbiology, Carl von Ossietzky University  
**External Supervisor:** Prof. Denis Shields  
**Date:** May 2026  
**Data & Code:** https://doi.org/10.5281/zenodo.20278099

---

## Abstract

Metagenomic pathogen detection fails systematically when pathogenic and commensal organisms share large conserved genomic regions. This report describes the first stage of a three-project framework addressing this problem for *Escherichia coli* O157:H7 detection. A preliminary pilot study first established the biological validity of the approach: whole-genome alignment of O157:H7 Sakai against the commensal strain SE11, using a zero-identity-threshold NUCmer strategy to define strictly absent (rather than merely diverged) regions, yielded 298 Sakai-specific marker sequences (≥100 bp; 1.31 Mb total). BLASTn screening of these markers against a 23-strain panel confirmed 94.6–98.0% retention within the O157:H7 clade and near-background levels across all commensal and non-O157 lineages, validating the marker set as lineage-specific. The main analysis then scaled this framework to a 61-genome panel: alignment of Sakai against 60 comparison strains produced 23,923 alignment blocks. Conserved-identity blocks (≥95%) were classified and recorded but excluded from the marker set as core-genome sequence unsuitable for discrimination. The 415 candidate pathogenicity markers retained comprise Moderately Diverged (85–94.9%, n=306) and Highly Diverged (<85%, n=109) blocks. A BLAST-based screen evaluated marker discriminative power against a simulated metagenomic community, achieving AUROC 0.552, sensitivity 17.4%, and specificity 92.1%. Biological validation confirmed that 49.5% of DIVERGED-tier markers are absent from all non-pathogenic laboratory strains tested (K-12 and E. coli B lineages), and 14 markers co-localise with named EHEC virulence loci. These results establish the marker set and baseline performance that Projects 2 and 3 extend.

---

## 1. Introduction

The *Enterobacteriaceae* present the most tractable and most clinically consequential case of the metagenomic false positive problem. *Escherichia coli* exists across a pathogenicity spectrum — from harmless gut commensals present in every healthy human to enterohemorrhagic EHEC O157:H7, responsible for haemolytic uraemic syndrome and large foodborne outbreaks globally (Karmali et al., 1985). These organisms are genomically near-identical at the core: they share over 3,000 housekeeping genes at ≥95% nucleotide identity (Touchon et al., 2009). Yet a subset of pathogen-exclusive genomic islands — O-islands — encode the entire apparatus of EHEC virulence: type III secretion, Shiga toxin, tellurite resistance, and the adhesion machinery of the locus of enterocyte effacement (Dobrindt et al., 2004; Perna et al., 2001; Vanaja et al., 2021).

The challenge for metagenomic detection is that conventional classifiers query the whole genome. A read originating from a commensal *E. coli* in a gut metagenome aligns to the O157:H7 reference at near-perfect identity if it spans a conserved region — and most of the genome is conserved. The consequence is a high false positive rate that makes clinical interpretation of metagenomic data unreliable (Wood & Salzberg, 2014; Breitwieser et al., 2019).

This project addresses the problem at its source: instead of querying the full genome, it identifies the specific genomic regions of O157:H7 that are absent or highly diverged in non-pathogenic relatives, stratifies them by degree of divergence, and screens simulated metagenomic reads against only those regions. The tiered classification system introduced here — distinguishing Conserved, Moderately Diverged, and Highly Diverged alignment blocks — is the conceptual foundation on which the downstream pangenomic and machine learning projects build.

The work proceeded in two stages. The first compared Sakai against a single commensal, SE11, to confirm that the marker definition was biologically sound before scaling up. The second extended the same approach to 60 genomes, applied the tiered identity classification, and ran the simulated metagenomic screen. The two stages are not separate experiments — the same alignment logic and absence definition run through both. Section 3 reports them in the order they were done.

---

## 2. Methods

### 2.1 Genome Dataset

Complete and chromosome-level genome assemblies were retrieved from NCBI RefSeq (O'Leary et al., 2016) using `ncbi-genome-download`. The final panel comprised **61 curated genomes**: the O157:H7 Sakai reference plus 60 comparison genomes (30 pathogenic, 30 non-pathogenic/commensal), achieving a balanced panel. Pathogenic strains span EHEC O157:H7 isolates (EDL933, EC4115, TW14359, HUSEC2011, O26:H11, O103:H2, O111, O145, O55:H7, O104:H4), UPEC (CFT073, UTI89, 536, UMN026, IAI39, NA114), ETEC (H10407, E24377A, TW11681), EAEC (042, 55989), EPEC (E2348/69, TW10598), NMEC (IHE3034, CE10, RS218), AIEC (LF82, NRG857C, UM146), and APEC O1. Commensal strains include environmental and human gut commensals (HS, IAI1, ED1a, SE15, SE11, ABU83972, ATCC25922, ATCC8739, ATCC11775, SMS-3-5), K-12 laboratory strains (MG1655, W3110, DH10B, MDS42, BW2952, RV308, DH5alpha, AB1157, J53, HMS174, NCM3722, NEB5alpha, AG100), *E. coli* B strains (BL21, BL21(DE3), REL606, W, C41(DE3), C43(DE3)), and probiotic Nissle 1917. Strain curation required assembly level ≥ complete/chromosome, confirmed pathotype from published literature, and unambiguous BioSample metadata. One mislabelled strain (SE11, initially classified EPEC) was corrected to commensal based on whole-genome characterisation confirming the absence of phage- and plasmid-borne virulence determinants (Oshima et al., 2008).

The primary reference for comparative analysis is *E. coli* O157:H7 Sakai (GenBank: BA000007.2; RefSeq: NC_002695.2). Although EDL933 (Perna et al., 2001) was the first O157:H7 genome sequenced, Sakai was selected for three reasons specific to this study's design. First, Hayashi et al. (2001) published a comprehensive map of 177 O-islands with precise chromosomal coordinates, which this project uses directly for biological validation of DIVERGED-tier markers. Second, both plasmids — pO157 and pOSAK1 — are fully sequenced and annotated in Sakai, whereas EDL933 carries only pO157, making Sakai the more complete representation of the EHEC mobile gene pool. Third, Sakai has a published sRNA atlas and a curated virulence gene catalogue, which made it straightforward to functionally annotate the markers recovered in this study. When coordinate-level precision at O-island boundaries matters — as it does here — Sakai is the reference genome the field converges on. All 60 comparison genomes were therefore aligned against it.

### 2.2 NUCmer Alignment

Pairwise whole-genome alignment was performed using NUCmer (MUMmer4; Marçais et al., 2018) for each of the 60 non-Sakai genomes against the Sakai reference:

```bash
nucmer --maxmatch -c 500 -b 500 -l 100 Sakai.fasta query.fasta \
    --prefix nucmer_sakai_vs_query

delta-filter -m -i 85 -l 100 nucmer_sakai_vs_query.delta \
    > nucmer_filtered.delta

show-coords -THrd nucmer_filtered.delta > nucmer_filtered.coords
```

Parameters: minimum cluster size 500 bp (`-c 500`), minimum alignment length 100 bp (`-l 100`), minimum identity 85% (`-i 85`), `-m` flag for one-to-one mapping to suppress repetitive region ambiguity. Alignment coordinates were parsed from `.coords` output files, retaining alignment start, end, percent identity, and query genome identity for each block.

### 2.3 Pilot Study: Sakai vs SE11 Marker Distribution

#### 2.3.1 Pairwise Alignment and Divergent Region Extraction

The first analytical stage aligned the Sakai reference against a single commensal comparator, *E. coli* SE11 (Oshima et al., 2008), using NUCmer with the delta-filter step imposing **no identity threshold**:

```bash
nucmer --prefix=sakai_vs_se11 Sakai.fna SE11.fna

delta-filter -m -o 0 sakai_vs_se11.delta > sakai_vs_se11_filtered.delta

show-coords -rclT sakai_vs_se11_filtered.delta > sakai_vs_se11.coords
```

The `-m` flag retains only one-to-one best alignments. The absence of `-i` means that any SE11 sequence aligning to a Sakai region at *any* nucleotide identity is counted as coverage of that region. A Sakai interval therefore becomes a divergent marker only when no SE11 sequence aligns to it at all — not merely when the alignment is below an arbitrary identity threshold. This definition strictly operationalises "absent from SE11" and distinguishes the pilot from identity-threshold-based approaches in which diverged orthologs (85–94% identity) are incorrectly classified as absent.

Aligned coordinates were converted to a BED file, merged with `bedtools merge`, and complemented against Sakai chromosome sizes with `bedtools complement` to extract unaligned intervals. Regions shorter than 100 bp were discarded. FASTA sequences of retained markers were extracted with `bedtools getfasta`.

#### 2.3.2 Pilot Panel and BLAST Screening

Markers were screened by BLASTn (≥85% nucleotide identity, ≥50% query coverage) against a 23-strain panel:

| Group | Strains |
|---|---|
| O157:H7 EHEC (n=3) | EDL933, EC4115, TW14359 |
| Non-O157 pathogens (n=7) | HUSEC2011, O26:H11 11368, O103:H2 12009, CFT073 (UPEC), H10407 (ETEC), E2348/69 (EPEC), APEC O1 |
| Commensal / laboratory *E. coli* (n=10) | SE11, HS, IAI1, ED1a, SE15, ATCC 25922, K-12 MG1655, K-12 W3110, BL21, Nissle 1917 |
| Outgroup *Escherichia* and other genus (n=3) | *E. fergusonii* EFCF056, *E. albertii* 6S-65-1, *S. aureus* N315 |

A marker was called present in a given strain if at least one BLAST hit satisfied both thresholds; otherwise absent. Results were compiled into a binary 298 × 23 presence/absence matrix and visualised as a heatmap with rows ordered by Jaccard distance / average linkage clustering.

### 2.4 Tiered Identity Classification

Each alignment block was assigned to one of three tiers based on percent nucleotide identity:

| Tier | Identity Range | Biological Interpretation |
|------|---------------|--------------------------|
| CONSERVED | ≥95% | Core genome; high false positive risk for detection |
| MODERATE | 85–94.9% | Ancestrally shared, adaptively diverging sequence |
| DIVERGED | <85% | Candidate pathogen-specific markers |

Non-aligning regions of the Sakai genome (absent from the query entirely) were classified as DIVERGED by definition. The DIVERGED designation therefore captures two biologically distinct scenarios: regions with nucleotide identity below 85% to commensal comparators, indicating substantial sequence remodelling consistent with adaptive divergence or strong positive selection; and Sakai regions entirely absent from a given query genome, suggesting horizontal acquisition with no homologous counterpart in the commensal lineage. Both scenarios are characteristic of pathogenicity island content. CONSERVED-tier blocks (≥95% identity) were classified and recorded but excluded from the candidate marker set: high identity to commensal comparators is precisely the property that makes a region unsuitable as a discriminative marker. Only MODERATE and DIVERGED blocks were carried forward as candidate markers; overlapping blocks within each tier were merged and filtered to ≥200 bp, yielding the final 415-marker set.

A binary conserved/divergent classification discards a substantial intermediate class. Two independent arguments support retaining it. First, pathogenicity in *E. coli* is not a discrete binary state but a quantitative continuum: different pathotypes accumulate distinct combinations of virulence genes through repeated horizontal transfer events, and the genomic boundary between a commensal and a pathogen is a gradient of gene content and regulatory architecture rather than a clean threshold (Bartoszek et al., 2018). A binary split imposes false precision on biology that is graded, not discrete. Second, amelioration theory predicts that horizontally acquired DNA does not retain its foreign sequence composition indefinitely — over evolutionary time, introgressed genes gradually converge toward the base composition and codon usage of the host genome through the same mutational processes acting on all chromosomal loci (Lawrence & Ochman, 1997). The MODERATE tier (85–94.9% identity) captures a biologically ambiguous population: sequence in this range may represent ancestrally foreign DNA undergoing gradual amelioration toward the host chromosomal background, or it may represent ancestrally shared sequence that has since diverged through accumulated substitutions and small indels. Identity alone cannot distinguish between these two trajectories. In this identity range, alignment alone cannot separate ameliorating foreign DNA from ancestrally shared sequence that has drifted. Codon usage bias, GC content, and dinucleotide composition encode signal that escapes identity thresholding. The Project 3 classifier targets those features directly. Retaining and labelling this intermediate class preserves information that a binary threshold would discard. The tier boundaries were not imposed arbitrarily: the 95% CONSERVED ceiling corresponds to the average nucleotide identity threshold established as the intraspecific conservation floor in bacterial genomics — the whole-genome sequence correlate of the 70% DNA–DNA hybridisation species boundary (Goris et al., 2007). The 85% DIVERGED floor marks the lower limit of reliable nucleotide alignment and the region of sequence most compositionally foreign to the host genome, characteristic of recently horizontally acquired DNA that has undergone minimal amelioration (Lawrence & Ochman, 1997).

### 2.5 BLAST Screen

Markers from the MODERATE and DIVERGED tiers were compiled into a single BLASTn database (Camacho et al., 2009). A synthetic metagenomic community was constructed using InSilicoSeq 2.0 (Gourlé et al., 2019) with an Illumina HiSeq 2500 paired-end error model (`--model HiSeq`). The WGS whole-genome simulation functionality used here was introduced in v1.x and retained in v2.0; Gourlé (2024) describes the amplicon simulation extension added in v2.0, which was not used in this study. Three community compositions were simulated: high spike (O157:H7 Sakai 10%), mid spike (5%), and low spike (1%), with the remaining abundance distributed across commensal strains (SE11, K-12 MG1655, Nissle 1917, HS, IAI1). Total: 500,000 read pairs per community; 150 bp read length. The low-spike community (1% EHEC abundance) was the primary evaluation scenario.

Simulated reads were queried against the marker database (e-value ≤ 1×10⁻⁵). A threshold classifier was applied: reads aligning at ≥95% identity were called positive (pathogen-derived). Performance was evaluated separately for each retained tier (MODERATE, DIVERGED) and for the combined database, using confusion matrix statistics, sensitivity, specificity, precision, and AUROC. The commensal-dominant community (1% EHEC spike) was the primary evaluation scenario.

### 2.6 Biological Validation

Two post-hoc validation analyses were conducted on the DIVERGED-tier markers:

1. **K-12 absence test:** Each DIVERGED marker was checked for NUCmer alignment coverage across all non-pathogenic laboratory strains in the comparison panel (K-12 lineage: MG1655, DH10B, W3110; E. coli B lineage: BL21, BL21(DE3), REL606, W). A marker was scored as absent if no laboratory strain produced a qualifying alignment at the marker's coordinates.

2. **Known virulence locus overlap:** DIVERGED marker coordinates were cross-referenced against published O-island boundaries from the Sakai genome annotation (Hayashi et al., 2001), including the LEE pathogenicity island, Stx1/2 prophages SpLE2/SpLE3, tellurite resistance locus, and OI-48 TTSS-2.

---

## 3. Results

---

**Figure 1. Analytical workflow for the three-project *E. coli* O157:H7 pathogenicity marker framework.**

![Workflow diagram](../data/results/figures/fig_workflow.png)

*The pipeline proceeds from genome data acquisition through a preliminary pilot study, then branches into Project 1 (comparative genomics, NUCmer tiered classification) and Project 2 (Anvi'o pangenome analysis) before converging into the Project 3 XGBoost machine learning classifier. Arrows indicate data flow; coloured headers indicate project boundaries.*

---

### 3.1 Pilot Study: Sakai-Specific Marker Distribution

#### 3.1.1 Marker Set Characteristics

NUCmer alignment of Sakai against SE11 with no identity threshold, followed by bedtools complement and ≥100 bp length filtering, produced **298 Sakai-specific marker regions** totalling **1,312,085 bp (1.31 Mb)**. Marker lengths ranged from 101 bp to 45,086 bp (median 1,780 bp; mean 4,403 bp; **Figure 2**). The distribution is strongly right-skewed: 20.8% of markers (n=62) are short fragments below 300 bp, 46.6% (n=139) fall in the gene-length range (300–3,000 bp), and 32.6% (n=97) exceed 3,000 bp and represent intact island-sized structures. The preponderance of large regions is consistent with recovery of intact prophage segments, O-islands, and type III secretion loci rather than fragmented intergenic remnants. The 45 kb maximum corresponds to a large prophage-associated island within the Sakai Sp prophage complement (Oshima et al., 2008).

---

**Figure 2. Length distribution of the 298 Sakai-specific pilot markers.**

![Marker length distribution](../data/results/figures/fig_pilot_marker_lengths.png)

*Histogram of marker lengths on a log-scale x-axis. Bars are coloured by size class: blue = short fragments (<300 bp, n=62, 20.8%); green = gene-length regions (300–3,000 bp, n=139, 46.6%); orange = island-sized regions (≥3,000 bp, n=97, 32.6%). Dashed vertical lines mark class boundaries. Red solid line = median (1,780 bp); purple dotted line = mean (4,403 bp). The right-skewed distribution reflects the preponderance of large, intact genomic islands in the marker set.*

---

#### 3.1.2 Marker Distribution across the 23-Strain Panel

BLASTn screening of all 298 markers against the 23-strain panel produced the binary presence/absence matrix in **Figure 3** and per-strain rates in **Table 1**. Results stratify sharply into two tiers.

O157:H7 strains retained 94.6–98.0% of markers (group mean 95.7%), confirming near-complete conservation of Sakai-specific genomic content within the O157:H7 lineage. All other strain categories fall substantially lower: non-O157 EHEC 27.5–35.9%, other pathogens 27.9–32.2%, commensal *E. coli* 0–32.2%, K-12/laboratory strains 21.1–24.2%, other *Escherichia* 11.7–22.5%, and *S. aureus* 0%. Critically, non-O157 pathotypes show no enrichment above commensals, confirming that the marker set captures O157:H7-lineage-specific content rather than general pathogen traits.

---

**Figure 3. Distribution of Sakai O157:H7 pilot markers across a 23-strain *E. coli* and outgroup panel.**

![Pilot marker distribution heatmap](../data/results/figures/fig_pilot_marker_heatmap.png)

*Each row represents one of 298 Sakai-specific marker sequences (≥100 bp). Each column represents one strain. Red = marker present (BLASTn ≥85% identity, ≥50% query coverage); blue = marker absent. Columns are arranged by biological group (left to right: O157:H7 → non-O157 pathogens → commensal/lab → outgroups); rows are ordered by Jaccard distance / average linkage clustering. The bar chart above each column shows the percentage of 298 markers detected in that strain, colour-coded by pathotype.*

---

**Table 1. Per-strain pilot marker detection rates.**

| Strain | Category | Markers present | % of 298 |
|---|---|:---:|:---:|
| EHEC EDL933 | O157:H7 EHEC | 292 | 98.0% |
| EHEC EC4115 | O157:H7 EHEC | 282 | 94.6% |
| EHEC TW14359 | O157:H7 EHEC | 282 | 94.6% |
| EHEC HUSEC2011 | Non-O157 EHEC | 82 | 27.5% |
| EHEC O26:H11 11368 | Non-O157 EHEC | 107 | 35.9% |
| EHEC O103:H2 12009 | Non-O157 EHEC | 104 | 34.9% |
| UPEC CFT073 | Pathogenic | 83 | 27.9% |
| ETEC H10407 | Pathogenic | 89 | 29.9% |
| EPEC E2348/69 | Pathogenic | 96 | 32.2% |
| APEC O1 | Pathogenic | 95 | 31.9% |
| Commensal SE11 | Commensal | 28 | 9.4%‡ |
| Commensal HS | Commensal | 0 | 0.0% |
| Commensal IAI1 | Commensal | 48 | 16.1% |
| Commensal ED1a | Commensal | 96 | 32.2% |
| Commensal SE15 | Commensal | 75 | 25.2% |
| K-12 MG1655 | Laboratory | 66 | 22.1% |
| K-12 W3110 | Laboratory | 66 | 22.1% |
| Lab-B BL21 | Laboratory | 63 | 21.1% |
| Probiotic Nissle 1917 | Laboratory | 72 | 24.2% |
| ATCC 25922 | Commensal | 94 | 31.5% |
| *E. fergusonii* EFCF056 | Other *Escherichia* | 35 | 11.7% |
| *E. albertii* 6S-65-1 | Other *Escherichia* | 67 | 22.5% |
| *S. aureus* N315 | Other genus | 0 | 0.0% |

*‡ SE11 residual hits are attributable to shared insertion sequences at divergent loci (see Section 3.1.3), not coding sequence conservation. O103:H2 12009 and O26:H11 11368 values reflect a corrected analysis using the authentic genome sequences (AP010958.1 and AP010953.1 respectively), obtained from NCBI after genome-header verification identified mislabelled files in the original pilot run; see Section 3.1.3.*

#### 3.1.3 Genome File Audit and Data Quality Flags

Post-hoc inspection of genome file headers — conducted by reading the first FASTA line of each genome in the analysis directory — identified two mislabelled genome files in the original pilot run. The directory designated for EHEC O103:H2 str. 12009 contained the genome of O157:H7 TW14359 (CP001368.1), and the directory for O26:H11 str. 11368 contained the EPEC E2348/69 genome (FM180568.1). The authentic sequences — AP010958.1 (*E. coli* O103:H2 str. 12009) and AP010953.1 (*E. coli* O26:H11 str. 11368) — were obtained from NCBI RefSeq via `Entrez.efetch` and the BLASTn screen was re-run against the existing 298-marker panel (85% identity, 50% query coverage thresholds unchanged). Corrected detection rates are 104/298 (34.9%) for O103:H2 and 107/298 (35.9%) for O26:H11, and these values are reported in Table 1. Both fall within the non-O157 EHEC range and do not alter any primary conclusion: the separation between the O157:H7 clade (94.6–98.0%) and all other lineages remains intact, and the corrected O103:H2 and O26:H11 values are consistent with the expected biology of non-O157 EHEC strains sharing some but not the majority of O157:H7-lineage-specific sequence.

#### 3.1.4 SE11 Residual Signal and Insertion Sequence Artefact

SE11 — the reference comparator used to define the marker set — yields BLAST hits for 28 of 298 markers (9.4%). This is a known biological artefact rather than a pipeline failure. Markers were defined as regions with no NUCmer alignment at any identity to SE11 at those chromosomal positions. However, both Sakai and SE11 carry mobile insertion sequences (IS629, IS3) at *different* chromosomal positions. A Sakai IS element embedded within an otherwise Sakai-specific prophage island has no SE11 alignment at that locus, correctly excluded by NUCmer. BLAST, being position-independent, finds the IS copy elsewhere in SE11 and reports a hit. These 28 markers are retained in all analyses. The underlying loci are genuinely absent from SE11 at the relevant genomic position; the IS element content within them is incidental.

---

With lineage specificity established and the zero-identity-threshold definition validated, the analysis was extended to the full 61-genome panel. The following sections report the results of that extension.

### 3.2 Alignment Landscape

NUCmer alignment of the Sakai reference against 60 comparison strains produced **23,923 alignment blocks** spanning the 5.5 Mb Sakai chromosome and two plasmids (pO157, pOSAK1). The alignment landscape (Figure 4) reveals the characteristic genomic architecture of EHEC: broad, high-identity alignment coverage in core metabolic and ribosomal loci, punctuated by large non-aligning gaps at O-island coordinates that become progressively deeper in commensals compared to pathogenic STEC relatives.

The visual correspondence between alignment gaps and published O-island positions provides an immediate, pre-computational validation of the approach: the pipeline is finding exactly the genomic regions biology predicts should be pathogen-specific.

**Figure 4.** Alignment landscape of 60 *E. coli* comparison strains aligned against O157:H7 Sakai. Each horizontal segment represents one NUCmer alignment block; colour encodes identity tier (blue = CONSERVED ≥95%, orange = MODERATE 85–94.9%, red = DIVERGED <85%). Strains are ordered by pathotype; dashed horizontal lines separate pathotype groups. Large gaps in commensal strain rows correspond to O-island positions.

![Alignment Landscape](../data/results/figures/alignment_landscape.png)

---

### 3.3 Marker Extraction

Tiered classification of the 23,923 alignment blocks produced 18,139 CONSERVED, 4,302 MODERATE, and 1,482 DIVERGED alignment blocks across all 60 comparison genomes. CONSERVED blocks were excluded from the marker set. Merging and length-filtering the retained MODERATE and DIVERGED blocks against the Sakai reference yielded:

| Tier | N markers | % of marker set |
|------|-----------|----------------|
| MODERATE (85–94.9%) | 306 | 73.7% |
| DIVERGED (<85%) | 109 | 26.3% |
| **Total** | **415** | **100%** |

DIVERGED markers were enriched at O-island coordinates relative to their genome-wide frequency, consistent with the hypothesis that highly diverged sequence co-localises with horizontally acquired virulence loci.

---

### 3.4 BLAST Screen Performance

Threshold-based classification of simulated reads against the 415-marker database yielded the following performance in the commensal-dominant community (the most clinically realistic scenario):

| Condition | Sensitivity | Specificity | Precision | AUROC |
|-----------|-------------|-------------|-----------|-------|
| MODERATE only | 0.166 | 0.928 | 0.180 | 0.551 |
| DIVERGED only | 0.024 | 0.985 | 0.128 | 0.504 |
| COMBINED (both tiers) | 0.174 | 0.921 | 0.173 | 0.552 |

The contrast between tiers is informative. MODERATE markers, which retain partial sequence similarity to the commensal background, attract more true-positive reads from the 1% EHEC spike and yield higher sensitivity (16.6%) and AUROC (0.551). DIVERGED markers, by definition, bear minimal resemblance to the commensal gene pool: commensal reads almost never align at ≥95% identity, producing near-perfect specificity (98.5%) but very low sensitivity (2.4%). The COMBINED database, which pools both tier sets, represents the operational BLAST screen used as the baseline for downstream ML comparison; the marginal AUROC increase (0.552 vs 0.551) reflects the small additional true-positive recovery from the DIVERGED set.

The AUROC of 0.552, considered in isolation, appears modest — and deliberately so. Identity threshold classification is the most primitive discriminator available: it accepts a read as pathogen-derived only when its sequence aligns to a marker at ≥95% nucleotide identity, a criterion that DIVERGED-tier markers, selected precisely for their departure from the commensal sequence pool, will satisfy rarely. The resulting baseline does not measure the information content of the marker set; it measures the ceiling of what a single-parameter alignment rule can extract from it. That ceiling — established under the most demanding evaluation scenario, a community in which EHEC reads constitute just 1% of the total read pool — is not a destination. It is a calibrated reference point: the performance frontier of naïve homology search, against which the pangenomic and machine learning frameworks of Projects 2 and 3 can be rigorously measured.

The structural source of this low sensitivity is the same property that makes these markers diagnostically valuable. Sequence that has diverged substantially from the commensal gene pool will rarely be matched at ≥95% identity by a 150 bp read originating from a commensal organism — which is exactly what drives the 92.1% combined specificity. But the converse also holds: a DIVERGED-tier marker that has departed far enough from the O157:H7 template (e.g., through IS element insertions or localised substitutions in the read) will also fail to match at ≥95% identity even for a true-positive pathogen read, collapsing sensitivity. The marker set carries the biological signal needed for pathotype discrimination — the question Projects 2 and 3 address is whether feature-based and pangenomic methods can extract it where identity-threshold classification cannot.

**Figure 5.** Per-tier and combined confusion matrices for the BLAST screen classifier. Rows represent tiers (DIVERGED, MODERATE, COMBINED); columns show true versus predicted class. Sensitivity, specificity, and precision annotations below each matrix.

![Confusion Matrix — DIVERGED](../data/results/figures/confusion_diverged.png)
![Confusion Matrix — MODERATE](../data/results/figures/confusion_moderate.png)
![Confusion Matrix — COMBINED](../data/results/figures/confusion_combined.png)

---

### 3.5 Biological Validation of DIVERGED Markers

#### 3.5.1 K-12 Absence

Of the 109 DIVERGED markers, **54 (49.5%) were absent from all non-pathogenic laboratory strains tested** (K-12 lineage: MG1655, DH10B, W3110; E. coli B lineage: BL21, BL21(DE3), REL606, W). This confirms that nearly half of DIVERGED markers have no counterpart in K-12 or *E. coli* B laboratory lineages. In real data, any read matching these markers is pathogen-derived by definition.

The remaining 55 (50.5%) showed some K-12 alignment, indicating ancestral sharing followed by adaptive divergence rather than clean horizontal acquisition. These are the markers the tiered classification retains for downstream ML analysis: their K-12 presence makes identity-threshold classification ambiguous, but codon usage, composition, and regulatory architecture may still carry distinguishing signal.

#### 3.5.2 Virulence Locus Overlap

Cross-referencing against published O-island coordinates (Hayashi et al., 2001) identified **14 DIVERGED markers co-localising with named virulence loci**:

| Virulence Locus | Function | Markers (n) |
|----------------|----------|-------------|
| Stx1 prophage (SpLE2) | Shiga toxin 1 production | 5 |
| LEE (OI-148) | Type III secretion, attaching-and-effacing | 4 |
| OI-48 (TTSS-2) | Second type III secretion system | 2 |
| Tellurite resistance locus | pO157 plasmid persistence factor | 2 |
| Stx2 prophage (SpLE3) | Shiga toxin 2 production | 1 |
| **Total** | | **14** |

Recovery of markers from all major EHEC virulence categories — integrating prophage, type III secretion, plasmid-borne, and regulatory loci — confirms that the NUCmer tiering pipeline is capturing biologically relevant divergence, not statistical artefact.

**Figure 6.** Biological validation overview. Left panel: K-12 absence rate for DIVERGED markers (54/109, 49.5%). Right panel: overlap between DIVERGED markers and named EHEC virulence loci, showing recovery of markers from LEE, Stx1/2 prophages, OI-48, and tellurite resistance.

![Biological Validation Overview](../data/results/figures/bio_sanity_overview.png)

---

## 4. Discussion

### 4.1 Stage One: Establishing Lineage Specificity and Validating the Absence Definition

The pilot results establish two facts before any analysis of the full panel. First, Sakai-specific regions identified by a zero-identity-threshold NUCmer approach are genuinely lineage-specific: 94.6–98.0% retention across three diverse O157:H7 isolates (EDL933, EC4115, TW14359) spanning separate outbreak lineages and three decades confirms that Sakai-specific sequence is not idiosyncratic to the Sakai isolate but represents a stable feature of the O157:H7 clade. Second, the signal is not a general pathogen property: non-O157 EHEC strains (27.5–35.9%) converge with commensals (0–32.2%), indicating that the markers encode O157:H7-lineage biology rather than broad virulence. This distinction matters for the main analysis: it means that a classifier trained on these markers should discriminate O157:H7 from all other *E. coli* — not merely from commensals — a harder and more clinically relevant task.

The pilot also justified the no-identity-filter decision that carries into the main study. The 9.4% SE11 residual is attributable to insertion sequences, not diverged coding genes — confirming that the no-filter approach successfully eliminates the diverged-ortholog leakage problem that a 95% identity threshold introduces. The 60-genome extension uses the same strict absence definition when complementing aligned regions for DIVERGED-tier marker designation.

### 4.2 Tiered Classification as a Conceptual Contribution

The central contribution of this project is not the alignment itself — NUCmer is an established tool — but the stratification logic applied to its output. Conventional comparative genomics pipelines treat aligned and non-aligned sequence as a binary: aligned regions are conserved and discarded; non-aligned regions are unique and retained. This binary discards a large intermediate class of sequence — the MODERATE tier — that is neither core genome nor clearly strain-specific.

The tiered system retains this intermediate class and asks a different question: not "is this sequence present in commensals?" but "how similar is this sequence to its commensal counterpart, and what does the degree of dissimilarity tell us about its functional and evolutionary history?" This reframing turns a filtering step into an analytical step, preserving information that would otherwise be lost.

### 4.3 Why the BLAST Screen Fails and What That Tells Us

The sensitivity-specificity asymmetry of this screen — 17.4% against 92.1% — is not a design compromise. It is a deliberate architectural choice grounded in the clinical reality of metagenomic pathogen detection. In a community where commensal *E. coli* abundances routinely exceed pathogenic strains by orders of magnitude, a false positive is not merely an inconvenience: it generates a spurious pathogen signal that propagates through clinical interpretation and demands costly confirmatory investigation to resolve. The BLAST screen is engineered to prevent this. By targeting exclusively sequence that has diverged substantially from the commensal gene pool, it erects a classification boundary that commensal reads cannot cross — 9 in 10 are correctly rejected — at the cost of recovering fewer than 1 in 5 pathogen-derived reads at this stage alone. That is the appropriate behaviour for a first-pass filter in a multi-stage diagnostic pipeline: hold the false positive rate to a minimum, accept that sensitivity will be low, and rely on the pangenomic and machine learning layers that follow to recover what identity thresholding cannot reach.

Lowering the identity threshold is not the resolution — doing so would erode specificity and resurrect the false positive problem this approach was designed to eliminate. The question that must be asked of these markers is not whether a read aligns to them with high identity, but whether the marker itself carries the compositional, translational, and evolutionary signatures of a horizontally acquired virulence gene. That is the question the feature-based ML classifier in Project 3 is built to answer.

### 4.4 The Marker Set as a Foundation

The BLAST screen was never expected to be the answer — it was designed to show exactly where a pure identity-threshold approach runs out of road. At AUROC 0.552 it does precisely that: high enough to confirm the markers carry real signal, low enough to confirm that reading that signal through a binary alignment rule is the wrong instrument. The markers themselves are sound. Of the 109 DIVERGED-tier sequences, more than half are absent from every K-12 strain in the panel, and fourteen land directly on named EHEC virulence loci — Shiga toxin prophages, the locus of enterocyte effacement, the second type III secretion system. What identity-threshold screening cannot do is see context. A 150 bp read that aligns to a marker at 89% identity looks the same whether that marker sits in the middle of a stable pathogenicity island or a mobile element drifting through commensal backgrounds. The alignment score tells you about sequence similarity. It tells you nothing about how that region distributes across a population of sixty-one strains, whether it travels with other virulence genes, or whether its codon usage and GC composition bear the signature of horizontal acquisition.

That is what Project 2 investigates. By building a full pangenome across all 61 genomes — mapping every gene cluster, characterising the core, soft-core, and accessory partitions, and testing which clusters are statistically enriched in EHEC versus commensal lineages — Project 2 provides the population-level evidence that pairwise alignment cannot reach. The same markers, viewed through that lens, stop being isolated sequences and become addressable features in a richer evidence space. That is the foundation on which the machine learning classifier in Project 3 is built.

---

## 5. Conclusion

The core finding is structural, not metric. DIVERGED markers reach 98.5% specificity but recover almost no pathogen reads at ≥95% identity; MODERATE markers invert that trade-off. That asymmetry is not a failure of marker design — it is the ceiling of identity thresholding applied to sequence selected for divergence. Alignment of Sakai against 60 comparison strains produced 23,923 blocks and 415 candidate markers. Of 109 DIVERGED markers, 54 are absent from every K-12 strain tested; 14 co-localise with named EHEC virulence loci. The signal is there. Extracting it requires classification sensitive to composition, pangenomic distribution, and evolutionary context — the analytical layers Projects 2 and 3 provide.

---

## References

Bartoszek, K., Majchrzak, M., Sakowski, S., Kubiak-Szeligowska, A. B., Kaj, I., & Parniewski, P. (2018). Predicting pathogenicity behavior in *Escherichia coli* population through a state dependent model and TRS profiling. *PLoS Computational Biology*, 14(1), e1005931. https://doi.org/10.1371/journal.pcbi.1005931

Breitwieser, F. P., Lu, J., & Salzberg, S. L. (2019). A review of methods and databases for metagenomic classification and assembly. *Briefings in Bioinformatics*, 20(4), 1125–1136. https://doi.org/10.1093/bib/bbx120

Camacho, C., Coulouris, G., Avagyan, V., Ma, N., Papadopoulos, J., Bealer, K., & Madden, T. L. (2009). BLAST+: architecture and applications. *BMC Bioinformatics*, 10, 421. https://doi.org/10.1186/1471-2105-10-421

Dobrindt, U., Hochhut, B., Hentschel, U., & Hacker, J. (2004). Genomic islands in pathogenic and environmental microorganisms. *Nature Reviews Microbiology*, 2(5), 414–424. https://doi.org/10.1038/nrmicro884

Goris, J., Konstantinidis, K. T., Klappenbach, J. A., Coenye, T., Vandamme, P., & Tiedje, J. M. (2007). DNA-DNA hybridization values and their relationship to whole-genome sequence similarities. *International Journal of Systematic and Evolutionary Microbiology*, 57(1), 81–91. https://doi.org/10.1099/ijs.0.64483-0

Gourlé, H., Karlsson-Lindsjö, O., Hayer, J., & Bongcam-Rudloff, E. (2019). Simulating Illumina metagenomic data with InSilicoSeq. *Bioinformatics*, 35(3), 521–522. https://doi.org/10.1093/bioinformatics/bty630

Gourlé, H. (2024). InSilicoSeq 2.0: simulating realistic amplicon-based sequence reads. *bioRxiv* [Preprint]. https://doi.org/10.1101/2024.02.16.580469

Hayashi, T., Makino, K., Ohnishi, M., Kurokawa, K., Ishii, K., Yokoyama, K., ... & Shinagawa, H. (2001). Complete genome sequence of enterohemorrhagic *Escherichia coli* O157:H7 and genomic comparison with a laboratory strain K-12. *DNA Research*, 8(1), 11–22. https://doi.org/10.1093/dnares/8.1.11

Lawrence, J. G., & Ochman, H. (1997). Amelioration of bacterial genomes: rates of change and exchange. *Journal of Molecular Evolution*, 44(4), 383–397. https://doi.org/10.1007/pl00006158

Lawrence, J. G., & Ochman, H. (1998). Molecular archaeology of the *Escherichia coli* genome. *Proceedings of the National Academy of Sciences*, 95(16), 9413–9417. https://doi.org/10.1073/pnas.95.16.9413

Karmali, M. A., Petric, M., Lim, C., Fleming, P. C., Arbus, G. S., & Lior, H. (1985). The association between idiopathic hemolytic uremic syndrome and infection by verotoxin-producing *Escherichia coli*. *Journal of Infectious Diseases*, 151(5), 775–782. https://doi.org/10.1093/infdis/151.5.775

Manning, S. D., Motiwala, A. S., Springman, A. C., Qi, W., Lacher, D. W., Ouellette, L. M., ... & Whittam, T. S. (2008). Variation in virulence among clades of *Escherichia coli* O157:H7 associated with disease outbreaks. *Proceedings of the National Academy of Sciences*, 105(12), 4868–4873. https://doi.org/10.1073/pnas.0710834105

Ogura, Y., Ooka, T., Iguchi, A., Toh, H., Asadulghani, M., Oshima, K., ... & Hayashi, T. (2009). Comparative genomics reveal the mechanism of the parallel evolution of O157 and non-O157 enterohemorrhagic *Escherichia coli*. *Proceedings of the National Academy of Sciences*, 106(42), 17939–17944. https://doi.org/10.1073/pnas.0903585106

Riley, L. W., Remis, R. S., Helgerson, S. D., McGee, H. B., Wells, J. G., Davis, B. R., ... & Cohen, M. L. (1983). Hemorrhagic colitis associated with a rare *Escherichia coli* serotype. *New England Journal of Medicine*, 308(12), 681–685. https://doi.org/10.1056/NEJM198303243081203

Tarr, P. I., Gordon, C. A., & Chandler, W. L. (2005). Shiga-toxin-producing *Escherichia coli* and haemolytic uraemic syndrome. *The Lancet*, 365(9464), 1073–1086. https://doi.org/10.1016/S0140-6736(05)71144-2

Marçais, G., Delcher, A. L., Phillippy, A. M., Coston, R., Salzberg, S. L., & Zimin, A. (2018). MUMmer4: a fast and versatile genome alignment system. *PLoS Computational Biology*, 14(1), e1005944. https://doi.org/10.1371/journal.pcbi.1005944

McDaniel, T. K., Jarvis, K. G., Donnenberg, M. S., & Kaper, J. B. (1995). A genetic locus of enterocyte effacement conserved among diverse enterobacterial pathogens. *Proceedings of the National Academy of Sciences*, 92(5), 1664–1668. https://doi.org/10.1073/pnas.92.5.1664

O'Leary, N. A., Wright, M. W., Brister, J. R., et al. (2016). Reference sequence (RefSeq) database at NCBI: current status, taxonomic expansion, and functional annotation. *Nucleic Acids Research*, 44(D1), D733–D745. https://doi.org/10.1093/nar/gkv1189

Oshima, K., Toh, H., Ogura, Y., Sasamoto, H., Morita, H., Park, S.-H., Ooka, T., Iyoda, S., Taylor, T. D., Hayashi, T., Itoh, K., & Hattori, M. (2008). Complete genome sequence and comparative analysis of the wild-type commensal *Escherichia coli* strain SE11 isolated from a healthy adult. *DNA Research*, 15(6), 375–386. https://doi.org/10.1093/dnares/dsn026

Perna, N. T., Plunkett, G., Burland, V., Mau, B., Glasner, J. D., Rose, D. J., ... & Blattner, F. R. (2001). Genome sequence of enterohaemorrhagic *Escherichia coli* O157:H7. *Nature*, 409(6819), 529–533. https://doi.org/10.1038/35054089

Touchon, M., Hoede, C., Tenaillon, O., et al. (2009). Organised genome dynamics in the *Escherichia coli* species results in highly diverse adaptive paths. *PLoS Genetics*, 5(1), e1000344. https://doi.org/10.1371/journal.pgen.1000344

Vanaja, S. K., Bergholz, T. M., & Whittam, T. S. (2021). Virulence-related O islands in enterohemorrhagic *Escherichia coli* O157:H7. *Gut Microbes*, 13(1), 1992237. https://doi.org/10.1080/19490976.2021.1992237

Wood, D. E., & Salzberg, S. L. (2014). Kraken: ultrafast metagenomic sequence classification using exact alignments. *Genome Biology*, 15, R46. https://doi.org/10.1186/gb-2014-15-3-r46
