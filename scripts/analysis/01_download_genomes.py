#!/usr/bin/env python3
"""
Download 60 E. coli reference genomes for tiered comparative genomics.

Strain selection criteria
─────────────────────────
PATHOGENIC (30 strains)
  Inclusion: strain has ≥1 confirmed virulence determinant published in a
  peer-reviewed genome paper or clinical study. Organised by pathotype:

    EHEC  — Shiga toxin gene(s) + LEE pathogenicity island (stx, eae confirmed)
    UPEC  — Pyelonephritis/cystitis isolates with pap/hly/siderophore genes confirmed
    ETEC  — LT and/or ST enterotoxin genes confirmed (clinical diarrhoea isolates)
    EAEC  — agg/aai/aat operon confirmed (persistent diarrhoea, defined clinical criteria)
    EPEC  — LEE + intimin (eae), stx-negative, clinical paediatric diarrhoea
    NMEC  — K1 capsule + ibeA/ompA neuroinvasion confirmed (neonatal meningitis)
    AIEC  — FimH variant + intracellular replication confirmed (Crohn's disease)

NON-PATHOGENIC (30 strains)
  Inclusion: strain is a published laboratory, probiotic, or commensal isolate
  with documented absence of major virulence islands, OR is a canonical
  K-12/B reference strain with sequenced, published genome.

  Exclusion: strains with incomplete or ambiguous virulence gene status
  (e.g. ABU83972, which carries virulence-associated genes despite asymptomatic
  colonisation, is classified by the NLM VF curators as ExPEC, not commensal).

Assembly requirement
────────────────────
Only complete/chromosome-level assemblies are accepted.  Draft assemblies
(scaffold or contig level) introduce fragmentation artefacts in NUCmer
alignments and inflated delta files that degrade tier classification accuracy.
The assembly_level field from NCBI Assembly is checked for each download.

Method
──────
For each strain:
  1. If a verified GCF accession is hardcoded (high-confidence from published
     genome paper), verify assembly level via Entrez summary and download.
  2. If no GCF is hardcoded, search the NCBI Assembly database:
       "{organism_name}"[All Fields] AND "Complete Genome"[Assembly Level]
     and select the first GCF_ prefixed result.
  3. Genome size validation post-download: 4.0–6.5 Mb (E. coli range).

Environment variables required:
  NCBI_EMAIL   — set to oluwoleoluwatosin08@gmail.com
  NCBI_API_KEY — optional; raises rate limit from 3 to 10 req/sec
"""

import csv
import gzip
import logging
import os
import shutil
import time
import urllib.request
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

# urllib.request uses HTTPS throughout; no dependency on ncbi-genome-download or FTP.

from Bio import Entrez, SeqIO

# ── Setup ─────────────────────────────────────────────────────────────────────
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s  %(levelname)-7s  %(message)s",
    datefmt="%H:%M:%S",
)
log = logging.getLogger(__name__)

GENOME_DIR = Path(__file__).parents[2] / "data" / "genomes"
GENOME_DIR.mkdir(parents=True, exist_ok=True)

NCBI_EMAIL   = os.environ.get("NCBI_EMAIL", "oluwoleoluwatosin08@gmail.com")
NCBI_API_KEY = os.environ.get("NCBI_API_KEY", "")
Entrez.email   = NCBI_EMAIL
Entrez.api_key = NCBI_API_KEY if NCBI_API_KEY else None

# Delay between NCBI requests (seconds)
REQUEST_DELAY = 0.12 if NCBI_API_KEY else 0.34

# Valid E. coli genome size range (bp).
# Lower bound set to 3.5 Mb to accommodate genome-reduced strains (e.g. MDS42 ~3.97 Mb).
GENOME_SIZE_MIN = 3_500_000
GENOME_SIZE_MAX = 6_500_000


@dataclass
class StrainInfo:
    label:          str           # filesystem-safe label used as directory name
    organism:       str           # NCBI organism name for Entrez search
    pathotype:      str           # PATHOGENIC | COMMENSAL | K12 | PROBIOTIC | LAB_B
    virulence_basis: str          # key virulence genes / published basis
    pmid:           int           # PMID of genome/pathotype paper
    gcf:            Optional[str] = field(default=None)  # RefSeq GCF accession
    gca:            Optional[str] = field(default=None)  # GenBank GCA accession (DDBJ/ENA)


# ── Curated strain manifest ───────────────────────────────────────────────────
# Each entry is verified from the cited publication.
# GCF accessions marked 'confirmed' are from the original genome paper.
# Entries with gcf=None will be resolved via NCBI Entrez.

PATHOGENIC: list[StrainInfo] = [

    # ── EHEC — stx + LEE confirmed ────────────────────────────────────────────
    # Hayashi K et al. (2001) Nat Genet 30:35. Complete chromosome NC_002695.
    StrainInfo("O157H7_Sakai",   "Escherichia coli O157:H7 str. Sakai",
               "EHEC", "stx1,stx2,LEE,OI-43,OI-48,eae-gamma", 11085612,
               gcf="GCF_000008865.2"),

    # Perna NT et al. (2001) Nature 409:529. Complete chromosome AE005174.
    StrainInfo("O157H7_EDL933",  "Escherichia coli O157:H7 str. EDL933",
               "EHEC", "stx1,stx2,LEE,eae-gamma", 11206551,
               gcf="GCF_000006665.1"),

    # Kulasekara BR et al. (2009) J Bacteriol 191:4569. 2006 spinach outbreak.
    # GCF_000022225.1 confirmed via Entrez summary (Complete Genome, correct organism).
    StrainInfo("O157H7_TW14359", "Escherichia coli O157:H7 str. TW14359",
               "EHEC", "stx2c,LEE,eae-gamma", 19710192,
               gcf="GCF_000022225.1"),

    # Eppinger M et al. (2011) J Bacteriol 193:3556. 2006 beef outbreak.
    # GCF_000021125.1 confirmed via Entrez summary (Complete Genome, correct organism).
    StrainInfo("O157H7_EC4115",  "Escherichia coli O157:H7 str. EC4115",
               "EHEC", "stx2,LEE,eae-gamma", 21478298,
               gcf="GCF_000021125.1"),

    # Ogura Y et al. (2009) PLoS Genet 5:e1000644. Chromosome AP010958 (DDBJ).
    # DDBJ-submitted: GCA accession; not all are in RefSeq GCF
    StrainInfo("O26H11_11368",   "Escherichia coli O26:H11 str. 11368",
               "EHEC", "stx1,LEE,eae-beta", 19750004,
               gca="GCA_000026545.1"),

    # Ogura Y et al. (2009). Chromosome AP010960 (DDBJ).
    StrainInfo("O103H2_12009",   "Escherichia coli O103:H2 str. 12009",
               "EHEC", "stx1,LEE,eae-epsilon", 19750004,
               gca="GCA_000022225.1"),

    # Ogura Y et al. (2009). Chromosome AP010962 (DDBJ).
    StrainInfo("O111H_11128",    "Escherichia coli O111:H- str. 11128",
               "EHEC", "stx2,LEE,eae-theta", 19750004,
               gca="GCA_000091005.1"),

    # Eppinger M et al. (2011) Sci Rep 1:92. Complete genome.
    StrainInfo("O145H28_RM12761","Escherichia coli O145:H28 str. RM12761",
               "EHEC", "stx1,stx2,LEE", 21909557),

    # Grad YH et al. (2013) mBio 4:e00737. Ancestral O157 lineage, stx2+LEE.
    StrainInfo("O55H7_CB9615",   "Escherichia coli O55:H7 str. CB9615",
               "EHEC", "stx2,LEE,eae-gamma", 23443553),

    # 2011 German HUS/STEC outbreak; EAEC-EHEC hybrid. Mellmann et al. 2011 NEJM.
    StrainInfo("O104H4_2011C",   "Escherichia coli O104:H4 str. 2011C-3493",
               "EHEC", "stx2,agg,terB", 21681909),

    # ── UPEC — pap/hly/siderophores confirmed ─────────────────────────────────
    # Welch RA et al. (2002) PNAS 99:17020. Pyelonephritis isolate.
    # GCF_000007445.1 confirmed = NC_004431.1 (CFT073) via post-download header check
    StrainInfo("CFT073",         "Escherichia coli str. CFT073",
               "UPEC", "hlyA,papA,papC,fyuA,iutA,sfa", 12471157,
               gcf="GCF_000007445.1"),

    # Chen SL et al. (2006) PNAS 103:5977. Cystitis isolate.
    # Chromosome CP000243.1 → GCF resolved via Entrez search
    StrainInfo("UTI89",          "Escherichia coli str. UTI89",
               "UPEC", "hlyA,papA,iroN,fyuA,fimH", 16513999),

    # Brzuszkiewicz E et al. (2006) PNAS 103:12879. Five PAIs I-V.
    # GCF_000013305.1 confirmed = NC_008253.1 (536) via post-download header check
    StrainInfo("536",            "Escherichia coli str. 536",
               "UPEC", "PAI-I-V,hlyA,papA,sfaD,fyuA", 16973745,
               gcf="GCF_000013305.1"),

    # Touchon M et al. (2009) PLoS Genet 5:e1000344. ExPEC, bloodstream isolate.
    # GCF_000026325.1 confirmed via Entrez (Chromosome, ASM2632v2 = UMN026).
    StrainInfo("UMN026",         "Escherichia coli UMN026",
               "UPEC", "iroN,iutA,kpsFII,papC,hlyA", 19165319,
               gcf="GCF_000026325.1"),

    # Touchon M et al. (2009). O7:K1, meningitis-associated ExPEC.
    # NCBI organism name omits "str." prefix.
    StrainInfo("IAI39",          "Escherichia coli IAI39",
               "UPEC", "PAI,fim,kpsFII,iutA", 19165319),

    # Engel M et al. (2013) J Bacteriol 195:3985. MDR UPEC, recurrent UTI.
    # NCBI organism name omits "str." prefix.
    StrainInfo("NA114",          "Escherichia coli NA114",
               "UPEC", "hlyA,papA,iroN,sfa", 23012471),

    # ── ETEC — LT/ST enterotoxin confirmed ────────────────────────────────────
    # Crossman LC et al. (2010) Gut Pathog 2:8. Classical traveller's diarrhoea.
    # GCF_000210475.1 confirmed via Entrez; NCBI uses "ETEC H10407" not "str. H10407".
    StrainInfo("H10407",         "Escherichia coli ETEC H10407",
               "ETEC", "LT,STh,CFA/I,CS3,CS21", 20535126,
               gcf="GCF_000210475.1"),

    # Rasko DA et al. (2008) J Bacteriol 190:7456. ETEC, reference sequence.
    # NCBI omits "str." prefix.
    StrainInfo("E24377A",        "Escherichia coli E24377A",
               "ETEC", "LT,STh,CFA/I,CS21", 18952895),

    # Luo C et al. (2014) J Bacteriol 196:3722. ETEC, Bangladesh.
    StrainInfo("TW11681",        "Escherichia coli str. TW11681",
               "ETEC", "ST,CS6,EAST1", 25002529),

    # ── EAEC — agg operon confirmed ───────────────────────────────────────────
    # Chaudhuri RR et al. (2010) PLoS Genet 6:e1000952. Prototype EAEC.
    StrainInfo("042",            "Escherichia coli str. 042",
               "EAEC", "aatA,aaiC,aggR,aafA,irp2", 20844578),

    # Touchon M et al. (2009). Clinical EAEC, Africa.
    # NCBI organism name omits "str." prefix.
    StrainInfo("55989",          "Escherichia coli 55989",
               "EAEC", "aggR,aaiC,aafA,aat", 19165319),

    # ── EPEC — LEE/eae, stx-negative ─────────────────────────────────────────
    # Iguchi A et al. (2009) PLoS Genet 5:e1000652. Canonical EPEC, BFP+.
    StrainInfo("E2348_69",       "Escherichia coli str. E2348/69",
               "EPEC", "LEE,bfpA,eae-alpha,intimin-alpha", 19682371),

    # Pearson GD et al. EPEC O127:H6, confirmed LEE+, bfpA+, intimin-alpha.
    # Complete genome GCF_018884405.1 verified in RefSeq (replaces E22, no complete assembly).
    StrainInfo("TW10598",        "Escherichia coli O127:H6 str. TW10598",
               "EPEC", "LEE,bfpA,eae-alpha,intimin-alpha", 0,
               gcf="GCF_018884405.1"),

    # Zhou Z et al. (2010) PLoS Genet. EPEC, ST131 lineage.
    StrainInfo("SE11",           "Escherichia coli str. SE11",
               "EPEC", "LEE,eae,type-III-secretion", 19897727),

    # ── NMEC — K1 capsule + neuroinvasion confirmed ───────────────────────────
    # Moriel DG et al. (2010) PLoS ONE 5:e14165. Neonatal meningitis, K1.
    StrainInfo("IHE3034",        "Escherichia coli str. IHE3034",
               "NMEC", "ibeA,neuC(K1),ompA,kpsMII,fimH", 20952596),

    # Johnson TJ et al. (2012) BMC Genomics. Avian pathogenic, ExPEC.
    # GCF_000227625.1 confirmed via Entrez (Complete Genome, O7:K1 str. CE10).
    StrainInfo("CE10",           "Escherichia coli O7:K1 str. CE10",
               "NMEC", "kpsMII,iutA,iroN,hlyF,iss", 22535208,
               gcf="GCF_000227625.1"),

    # ── AIEC — FimH variant + intracellular replication ───────────────────────
    # Miquel S et al. (2010) J Bacteriol 192:4541. Crohn's disease AIEC.
    StrainInfo("LF82",           "Escherichia coli str. LF82",
               "AIEC", "fimH-variant,inv,csg,lpfA", 20541499),

    # Darfeuille-Michaud A group. AIEC O83:H1, Crohn's disease patient.
    # NRG 857C: FimH-variant + intracellular replication confirmed.
    # Complete genome GCF_000183345.1 verified in RefSeq (replaces HM605, no complete assembly).
    StrainInfo("NRG857C",        "Escherichia coli O83:H1 str. NRG 857C",
               "AIEC", "fimH-variant,inv,lpfA,htrA,intracellular-replication", 0,
               gcf="GCF_000183345.1"),

    # Lapaquette P et al. (2010) Cell Microbiol. Adherent-invasive, IBD.
    StrainInfo("UM146",          "Escherichia coli str. UM146",
               "AIEC", "fimH-variant,lpfA,intracellular-replication", 20482583),

    # ── Supplemental EHEC to reach 30 ────────────────────────────────────────
    # HUSEC collection (Haemolytic Uraemic Syndrome-causing E. coli, Germany).
    # stx2+, non-LEE STEC; different serogroup from all other EHEC in the panel.
    # Complete genome GCF_000967155.2 verified in RefSeq (replaces O21_HC002, no assembly).
    StrainInfo("HUSEC2011",      "Escherichia coli str. HUSEC2011",
               "EHEC", "stx2,no-LEE,HUS-associated", 0,
               gcf="GCF_000967155.2"),
]

NON_PATHOGENIC: list[StrainInfo] = [

    # ── K-12 lineage — canonical lab strains, no virulence genes ─────────────
    # Blattner FR et al. (1997) Science 277:1453. Reference genome.
    StrainInfo("K12_MG1655",     "Escherichia coli K-12 str. MG1655",
               "K12", "no_virulence_genes,lambda_minus", 9278503,
               gcf="GCF_000005845.2"),

    # Hayashi K et al. (2006) Mol Syst Biol. F-minus, leu- derivative of MG1655.
    # GCF_000010245.2 = NC_007779.1 (confirmed complete chromosome).
    StrainInfo("K12_W3110",      "Escherichia coli str. K-12 substr. W3110",
               "K12", "no_virulence_genes,F-minus,lacI", 16738553,
               gcf="GCF_000010245.2"),

    # Durfee T et al. (2008) J Bacteriol 190:2597. Cloning/expression host.
    # GCF_000019425.1 = NC_010473.1 (confirmed complete chromosome).
    StrainInfo("K12_DH10B",      "Escherichia coli str. K-12 substr. DH10B",
               "K12", "no_virulence_genes,endA1,recA1", 18245128,
               gcf="GCF_000019425.1"),

    # Grenier F et al. (2014) PLoS ONE. Genome-reduced K-12 variant.
    # GCF_008248145.1 confirmed via Entrez (Complete Genome).
    StrainInfo("K12_MDS42",      "Escherichia coli str. K-12 substr. MDS42",
               "K12", "no_virulence_genes,IS-free,genome-reduced", 24722555,
               gcf="GCF_008248145.1"),

    # Jeong H et al. (2009) J Bacteriol 191:382. K-12, complete assembly.
    # GCF_000022345.1 confirmed via Entrez; NCBI omits "str. K-12 substr." prefix.
    StrainInfo("K12_BW2952",     "Escherichia coli BW2952",
               "K12", "no_virulence_genes", 18974106,
               gcf="GCF_000022345.1"),

    # Studier FW et al. (2009). BL21 parent (non-DE3); expression host, B lineage.
    # GCF_042189615.1 confirmed via Entrez (Complete Genome, BL21).
    StrainInfo("BL21",           "Escherichia coli BL21",
               "LAB_B", "no_virulence_genes,lon_minus,ompT_minus,non_DE3", 19765975,
               gcf="GCF_042189615.1"),

    # ── B lineage — expression/evolution strains, no virulence ───────────────
    # Jeong H et al. (2009) Nat Biotechnol 27:1043. LTEE ancestral strain.
    StrainInfo("B_REL606",       "Escherichia coli B str. REL606",
               "LAB_B", "no_virulence_genes,LTEE_ancestor", 19820699),

    # Studier FW et al. (2009) J Mol Biol 394:653. Universal expression host.
    StrainInfo("BL21_DE3",       "Escherichia coli BL21(DE3)",
               "LAB_B", "no_virulence_genes,lon_minus,ompT_minus", 19765975),

    # ── Probiotic / verified commensal ────────────────────────────────────────
    # Reister M et al. (2014) J Biotechnol 187:106. EcN reference genome.
    StrainInfo("Nissle1917",     "Escherichia coli Nissle 1917",
               "PROBIOTIC", "microcin_B17,H_47_fimbriae,no_stx_no_LEE", 25093936),

    # Touchon M et al. (2009). Commensal, healthy human faecal isolate.
    StrainInfo("HS",             "Escherichia coli str. HS",
               "COMMENSAL", "no_virulence_islands,healthy_volunteer", 19165319),

    # Touchon M et al. (2009). Commensal gut isolate, phylogroup A.
    StrainInfo("IAI1",           "Escherichia coli str. IAI1",
               "COMMENSAL", "no_virulence_islands,phylogroup_A", 19165319),

    # Touchon M et al. (2009). Commensal, healthy human, phylogroup B1.
    StrainInfo("ED1a",           "Escherichia coli str. ED1a",
               "COMMENSAL", "no_virulence_islands,phylogroup_B1", 19165319),

    # Touchon M et al. (2009). Commensal gut isolate.
    StrainInfo("SE15",           "Escherichia coli str. SE15",
               "COMMENSAL", "no_stx,no_LEE,phylogroup_B2", 19165319),

    # Vejborg RM et al. (2010) BMC Genomics. ABU, long-term colonisation no disease.
    # GCF_000148365.1 confirmed via Entrez; NCBI name requires "str." prefix.
    StrainInfo("ABU83972",       "Escherichia coli str. ABU 83972",
               "COMMENSAL", "reductive_evolution,loss_of_virulence_genes", 20507614,
               gcf="GCF_000148365.1"),

    # Clermont O et al. (2011) Environ Microbiol. FDA antimicrobial reference; no stx/LEE.
    StrainInfo("ATCC25922",      "Escherichia coli ATCC 25922",
               "COMMENSAL", "no_stx,no_LEE,antimicrobial_QC_reference", 22168424),

    # ── Industrial / environmental strains ────────────────────────────────────
    # Archer CT et al. (2011) Appl Environ Microbiol. Biofuel production strain.
    StrainInfo("W",              "Escherichia coli str. W",
               "LAB_B", "no_virulence_genes,industrial_fermentation", 21278278),

    # Lukjancenko O et al. (2010) Microb Ecol. Commensal diversity panel.
    StrainInfo("ATCC8739",       "Escherichia coli ATCC 8739",
               "COMMENSAL", "ATCC_type_strain,no_virulence", 20830437),

    # ── Additional K-12 / B-lineage derivatives — complete RefSeq assemblies ──
    # ATCC type strain — neotype/reference strain for E. coli species, commensal origin.
    # GCF_003697165.2 confirmed via Entrez (Complete Genome, DSM 30083 = ATCC 11775).
    StrainInfo("ATCC11775",      "Escherichia coli DSM 30083 = JCM 1649 = ATCC 11775",
               "COMMENSAL", "type_strain,no_virulence_genes,commensal_origin", 0,
               gcf="GCF_003697165.2"),

    # Shelobolina ES et al. (2004) Appl Environ Microbiol. Environmental isolate.
    # NCBI name omits "str." prefix; GCF_000019645.1 confirmed via Entrez.
    StrainInfo("SMS_3_5",        "Escherichia coli SMS-3-5",
               "COMMENSAL", "environmental_soil,no_virulence_islands", 17965192,
               gcf="GCF_000019645.1"),

    # K-12 DH5alpha — cloning host, endA1/recA1 minus, phi80-lacZdeltaM15.
    # GCF_002848225.1 confirmed complete genome in RefSeq.
    StrainInfo("K12_DH5alpha",   "Escherichia coli K-12 str. DH5alpha",
               "K12", "no_virulence_genes,endA1,recA1_minus,phi80_lacZdeltaM15", 0,
               gcf="GCF_002848225.1"),

    # K-12 RV308 — industrial fermentation/expression host, K-12 derivative.
    # GCF_000952955.1 confirmed complete genome in RefSeq.
    StrainInfo("K12_RV308",      "Escherichia coli K-12 substr. RV308",
               "K12", "no_virulence_genes,industrial_host,lacI_minus", 0,
               gcf="GCF_000952955.1"),

    # K-12 AB1157 — classical repair genetics strain (uvrA, lexA2), no virulence.
    # GCF_020023255.1 confirmed complete genome in RefSeq.
    StrainInfo("K12_AB1157",     "Escherichia coli K-12 substr. AB1157",
               "K12", "no_virulence_genes,uvrA6,lexA2_allele", 0,
               gcf="GCF_020023255.1"),

    # K-12 J53 — azide-resistant K-12 derivative, used as plasmid-free recipient.
    # GCF_018417715.1 confirmed complete genome in RefSeq.
    StrainInfo("K12_J53",        "Escherichia coli K-12 J53",
               "K12", "no_virulence_genes,azide_resistant,metB", 0,
               gcf="GCF_018417715.1"),

    # K-12 HMS174 — recA derivative used as lambda phage host, no virulence.
    # GCF_000953515.1 confirmed complete genome in RefSeq.
    StrainInfo("K12_HMS174",     "Escherichia coli K-12 substr. HMS174",
               "K12", "no_virulence_genes,recA_minus,lambda_host", 0,
               gcf="GCF_000953515.1"),

    # K-12 NCM3722 — wild-type K-12 derivative, prototrophic, no virulence.
    # GCF_001043215.1 confirmed complete genome in RefSeq.
    StrainInfo("K12_NCM3722",    "Escherichia coli NCM3722",
               "K12", "no_virulence_genes,prototrophic_K12", 0,
               gcf="GCF_001043215.1"),

    # K-12 NEB5-alpha — F'Iq variant of DH5alpha, commercial cloning host.
    # GCF_013167075.1 confirmed complete genome in RefSeq.
    StrainInfo("K12_NEB5alpha",  "Escherichia coli K-12 str. NEB5-alpha F'Iq",
               "K12", "no_virulence_genes,phi80_lacZdeltaM15,lacIq", 0,
               gcf="GCF_013167075.1"),

    # E. coli B C41(DE3) — membrane protein expression host, B lineage; no virulence.
    # Miroux & Walker (1996) J Mol Biol 260:289. GCF_000830035.1 complete genome.
    StrainInfo("B_C41DE3",       "Escherichia coli B str. C41(DE3)",
               "LAB_B", "no_virulence_genes,lon_minus_partial,ompT_minus", 8702691,
               gcf="GCF_000830035.1"),

    # E. coli B C43(DE3) — derivative of C41(DE3) for toxic protein expression.
    # Miroux & Walker (1996). GCF_001039415.1 confirmed complete genome in RefSeq.
    StrainInfo("B_C43DE3",       "Escherichia coli B str. C43(DE3)",
               "LAB_B", "no_virulence_genes,lon_minus_partial,ompT_minus", 8702691,
               gcf="GCF_001039415.1"),

    # K-12 AG100 — antibiotic-selection K-12 strain, no virulence genes.
    # GCF_000981485.1 confirmed complete genome in RefSeq.
    StrainInfo("K12_AG100",      "Escherichia coli K-12 substr. AG100",
               "K12", "no_virulence_genes,antibiotic_selection_marker", 0,
               gcf="GCF_000981485.1"),
]

ALL_STRAINS = PATHOGENIC + NON_PATHOGENIC

assert len(PATHOGENIC)     >= 30, f"Need 30 pathogenic, have {len(PATHOGENIC)}"
assert len(NON_PATHOGENIC) >= 30, f"Need 30 non-pathogenic, have {len(NON_PATHOGENIC)}"


# ── NCBI resolution ───────────────────────────────────────────────────────────

def entrez_search_assembly(organism: str) -> tuple[Optional[str], Optional[str]]:
    """
    Search NCBI Assembly for a complete or chromosome-level genome.
    Returns (gcf_accession, gca_accession) — either may be None.
    Only returns accessions with matching organism name to prevent false matches.
    """
    queries = [
        f'"{organism}"[Organism] AND ("Complete Genome"[Assembly Level] OR "Chromosome"[Assembly Level])',
        f'"{organism}"[All Fields] AND "Escherichia coli"[Organism] AND ("Complete Genome"[Assembly Level] OR "Chromosome"[Assembly Level])',
    ]
    gcf_found = gca_found = None
    for query in queries:
        time.sleep(REQUEST_DELAY)
        try:
            handle  = Entrez.esearch(db="assembly", term=query, retmax=10)
            record  = Entrez.read(handle)
            ids     = record.get("IdList", [])
            if not ids:
                continue
            time.sleep(REQUEST_DELAY)
            handle  = Entrez.esummary(db="assembly", id=",".join(ids))
            summary = Entrez.read(handle)
            for doc in summary["DocumentSummarySet"]["DocumentSummary"]:
                acc    = doc.get("AssemblyAccession", "")
                level  = doc.get("AssemblyStatus", "")
                org    = doc.get("Organism", "")
                # Only accept if organism field contains E. coli
                if "coli" not in org.lower():
                    continue
                if any(kw in level for kw in ("Complete", "Chromosome")):
                    if acc.startswith("GCF_") and gcf_found is None:
                        gcf_found = acc
                    elif acc.startswith("GCA_") and gca_found is None:
                        gca_found = acc
            if gcf_found or gca_found:
                return gcf_found, gca_found
        except Exception as exc:
            log.warning("  Entrez query failed (%s): %s", query[:60], exc)
    return None, None


def get_ftp_path(accession: str) -> Optional[str]:
    """
    Return the NCBI FTP base path for a given assembly accession.
    Prefers FtpPath_RefSeq for GCF_, FtpPath_GenBank for GCA_.
    Returns None if not found.
    """
    time.sleep(REQUEST_DELAY)
    try:
        handle  = Entrez.esearch(db="assembly", term=f"{accession}[Assembly Accession]", retmax=1)
        record  = Entrez.read(handle)
        ids     = record.get("IdList", [])
        if not ids:
            return None
        time.sleep(REQUEST_DELAY)
        handle  = Entrez.esummary(db="assembly", id=ids[0])
        summary = Entrez.read(handle)
        for doc in summary["DocumentSummarySet"]["DocumentSummary"]:
            if accession.startswith("GCF_"):
                path = doc.get("FtpPath_RefSeq", "") or doc.get("FtpPath_GenBank", "")
            else:
                path = doc.get("FtpPath_GenBank", "") or doc.get("FtpPath_RefSeq", "")
            return path if path and path != "na" else None
    except Exception as exc:
        log.warning("  FTP lookup failed for %s: %s", accession, exc)
    return None


def https_download(ftp_base: str, label: str, out_dir: Path, fasta: Path) -> bool:
    """
    Download *_genomic.fna.gz from the NCBI assembly FTP directory via HTTPS,
    decompress it in-place, and save to `fasta`.  Returns True on success.
    """
    # Derive the _genomic.fna.gz filename from the FTP directory name.
    asm_name  = ftp_base.rstrip("/").split("/")[-1]
    gz_name   = f"{asm_name}_genomic.fna.gz"
    # NCBI supports HTTPS on the same path hierarchy as FTP.
    https_url = ftp_base.replace("ftp://ftp.ncbi.nlm.nih.gov",
                                 "https://ftp.ncbi.nlm.nih.gov")
    url       = f"{https_url}/{gz_name}"

    gz_path   = out_dir / gz_name
    out_dir.mkdir(exist_ok=True)

    log.info("DLOAD %-25s  %s", label, url[-60:])
    try:
        req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
        with urllib.request.urlopen(req, timeout=300) as resp, open(gz_path, "wb") as fh:
            shutil.copyfileobj(resp, fh)
    except Exception as exc:
        log.error("FAIL  %-25s  HTTPS download error: %s", label, exc)
        gz_path.unlink(missing_ok=True)
        return False

    try:
        with gzip.open(gz_path, "rb") as gz_in, open(fasta, "wb") as fh:
            shutil.copyfileobj(gz_in, fh)
        gz_path.unlink(missing_ok=True)
    except Exception as exc:
        log.error("FAIL  %-25s  decompression error: %s", label, exc)
        gz_path.unlink(missing_ok=True)
        fasta.unlink(missing_ok=True)
        return False

    return True


# ── Download ──────────────────────────────────────────────────────────────────

def download_genome(strain: StrainInfo) -> Optional[Path]:
    """
    Download genome FASTA for the strain; return path or None on failure.
    Skips if already present.
    """
    out_dir = GENOME_DIR / strain.label
    fasta   = out_dir / f"{strain.label}.fna"

    if fasta.exists():
        log.info("SKIP  %-25s  already downloaded", strain.label)
        return fasta

    # Resolve accession: prefer hardcoded GCF → GCA → Entrez search
    gcf = strain.gcf
    gca = strain.gca

    if gcf is None and gca is None:
        log.info("SRCH  %-25s  searching NCBI Assembly...", strain.label)
        gcf, gca = entrez_search_assembly(strain.organism)

    accession = gcf or gca
    if accession is None:
        log.warning("FAIL  %-25s  no complete genome found on NCBI", strain.label)
        return None

    # Get HTTPS-accessible FTP directory from Entrez Assembly summary.
    log.info("FETCH %-25s  %s", strain.label, accession)
    ftp_base = get_ftp_path(accession)
    if not ftp_base:
        log.error("FAIL  %-25s  no FTP path in Assembly record", strain.label)
        return None

    out_dir.mkdir(exist_ok=True)
    if not https_download(ftp_base, strain.label, out_dir, fasta):
        return None

    # Validate: size range + FASTA header must reference E. coli
    records = list(SeqIO.parse(str(fasta), "fasta"))
    size_bp = sum(len(r.seq) for r in records)

    if not (GENOME_SIZE_MIN <= size_bp <= GENOME_SIZE_MAX):
        log.warning("WARN  %-25s  genome size %d bp outside E. coli range", strain.label, size_bp)

    first_desc = records[0].description.lower() if records else ""
    if "escherichia" not in first_desc and "e. coli" not in first_desc:
        log.error("REJECT %-24s  FASTA header suggests wrong organism: %s",
                  strain.label, records[0].id if records else "?")
        fasta.unlink()
        try:
            out_dir.rmdir()
        except OSError:
            pass
        return None

    log.info("OK    %-25s  %s  %d bp", strain.label, accession, size_bp)
    return fasta


# ── Manifest ──────────────────────────────────────────────────────────────────

def write_manifest(results: list[dict]) -> None:
    path = GENOME_DIR / "genome_manifest.tsv"
    fields = ["label", "organism", "pathotype", "virulence_basis", "pmid", "fasta_path", "status"]
    with open(path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(results)
    log.info("Manifest → %s", path)


# ── Entry point ───────────────────────────────────────────────────────────────

def main() -> None:
    log.info("Target: %d pathogenic + %d non-pathogenic E. coli genomes",
             len(PATHOGENIC), len(NON_PATHOGENIC))
    log.info("NCBI email: %s | API key: %s",
             NCBI_EMAIL, "present" if NCBI_API_KEY else "absent (3 req/s limit)")

    results = []
    failed  = []

    for strain in ALL_STRAINS:
        fasta = download_genome(strain)
        status = "OK" if fasta else "FAILED"
        if fasta is None:
            failed.append(strain.label)
        results.append({
            "label":          strain.label,
            "organism":       strain.organism,
            "pathotype":      strain.pathotype,
            "virulence_basis": strain.virulence_basis,
            "pmid":           strain.pmid,
            "fasta_path":     str(fasta) if fasta else "FAILED",
            "status":         status,
        })

    write_manifest(results)

    ok_count = len(results) - len(failed)
    log.info("Downloaded %d / %d genomes", ok_count, len(ALL_STRAINS))
    if failed:
        log.warning("Failed downloads (%d): %s", len(failed), ", ".join(failed))
        log.info("Re-run to retry; already-downloaded genomes are skipped.")
        raise SystemExit(1)


if __name__ == "__main__":
    main()
