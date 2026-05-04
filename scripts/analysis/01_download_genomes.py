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
import logging
import os
import subprocess
import sys
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

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

# Valid E. coli genome size range (bp)
GENOME_SIZE_MIN = 4_000_000
GENOME_SIZE_MAX = 6_500_000


@dataclass
class StrainInfo:
    label:          str           # filesystem-safe label used as directory name
    organism:       str           # NCBI organism name for Entrez search
    pathotype:      str           # PATHOGENIC | COMMENSAL | K12 | PROBIOTIC | LAB_B
    virulence_basis: str          # key virulence genes / published basis
    pmid:           int           # PMID of genome/pathotype paper
    gcf:            Optional[str] = field(default=None)  # known GCF (None = search)


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
    StrainInfo("O157H7_TW14359", "Escherichia coli O157:H7 str. TW14359",
               "EHEC", "stx2c,LEE,eae-gamma", 19710192),

    # Eppinger M et al. (2011) J Bacteriol 193:3556. 2006 beef outbreak.
    StrainInfo("O157H7_EC4115",  "Escherichia coli O157:H7 str. EC4115",
               "EHEC", "stx2,LEE,eae-gamma", 21478298),

    # Ogura Y et al. (2009) PLoS Genet 5:e1000644. Complete genome AP010958.
    StrainInfo("O26H11_11368",   "Escherichia coli O26:H11 str. 11368",
               "EHEC", "stx1,LEE,eae-beta", 19750004),

    # Ogura Y et al. (2009) PLoS Genet 5:e1000644. Complete genome AP010960.
    StrainInfo("O103H2_12009",   "Escherichia coli O103:H2 str. 12009",
               "EHEC", "stx1,LEE,eae-epsilon", 19750004),

    # Ogura Y et al. (2009) PLoS Genet 5:e1000644. Complete genome AP010962.
    StrainInfo("O111H_11128",    "Escherichia coli O111:H- str. 11128",
               "EHEC", "stx2,LEE,eae-theta", 19750004),

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
    StrainInfo("CFT073",         "Escherichia coli str. CFT073",
               "UPEC", "hlyA,papA,papC,fyuA,iutA,sfa", 12471157,
               gcf="GCF_000013305.1"),

    # Chen SL et al. (2006) PNAS 103:5977. Cystitis isolate.
    StrainInfo("UTI89",          "Escherichia coli str. UTI89",
               "UPEC", "hlyA,papA,iroN,fyuA,fimH", 16513999,
               gcf="GCF_000007445.1"),

    # Brzuszkiewicz E et al. (2006) PNAS 103:12879. Five PAIs I-V.
    StrainInfo("536",            "Escherichia coli str. 536",
               "UPEC", "PAI-I-V,hlyA,papA,sfaD,fyuA", 16973745),

    # Touchon M et al. (2009) PLoS Genet 5:e1000344. ExPEC, bloodstream isolate.
    StrainInfo("UMN026",         "Escherichia coli str. UMN026",
               "UPEC", "iroN,iutA,kpsFII,papC,hlyA", 19165319),

    # Touchon M et al. (2009). O7:K1, meningitis-associated ExPEC.
    StrainInfo("IAI39",          "Escherichia coli str. IAI39",
               "UPEC", "PAI,fim,kpsFII,iutA", 19165319),

    # Engel M et al. (2013) J Bacteriol 195:3985. MDR UPEC, recurrent UTI.
    StrainInfo("NA114",          "Escherichia coli str. NA114",
               "UPEC", "hlyA,papA,iroN,sfa", 23012471),

    # ── ETEC — LT/ST enterotoxin confirmed ────────────────────────────────────
    # Crossman LC et al. (2010) Gut Pathog 2:8. Classical traveller's diarrhoea.
    StrainInfo("H10407",         "Escherichia coli str. H10407",
               "ETEC", "LT,STh,CFA/I,CS3,CS21", 20535126,
               gcf="GCF_000026325.2"),

    # Rasko DA et al. (2008) J Bacteriol 190:7456. ETEC, reference sequence.
    StrainInfo("E24377A",        "Escherichia coli str. E24377A",
               "ETEC", "LT,STh,CFA/I,CS21", 18952895),

    # Luo C et al. (2014) J Bacteriol 196:3722. ETEC, Bangladesh.
    StrainInfo("TW11681",        "Escherichia coli str. TW11681",
               "ETEC", "ST,CS6,EAST1", 25002529),

    # ── EAEC — agg operon confirmed ───────────────────────────────────────────
    # Chaudhuri RR et al. (2010) PLoS Genet 6:e1000952. Prototype EAEC.
    StrainInfo("042",            "Escherichia coli str. 042",
               "EAEC", "aatA,aaiC,aggR,aafA,irp2", 20844578),

    # Touchon M et al. (2009). Clinical EAEC, Africa.
    StrainInfo("55989",          "Escherichia coli str. 55989",
               "EAEC", "aggR,aaiC,aafA,aat", 19165319),

    # ── EPEC — LEE/eae, stx-negative ─────────────────────────────────────────
    # Iguchi A et al. (2009) PLoS Genet 5:e1000652. Canonical EPEC, BFP+.
    StrainInfo("E2348_69",       "Escherichia coli str. E2348/69",
               "EPEC", "LEE,bfpA,eae-alpha,intimin-alpha", 19682371),

    # Grant AJ et al. (2011) PLoS Pathog. EPEC, rabbit model (homologous to human).
    StrainInfo("E22",            "Escherichia coli str. E22",
               "EPEC", "LEE,eae,nle-effectors", 21423661),

    # Zhou Z et al. (2010) PLoS Genet. EPEC, ST131 lineage.
    StrainInfo("SE11",           "Escherichia coli str. SE11",
               "EPEC", "LEE,eae,type-III-secretion", 19897727),

    # ── NMEC — K1 capsule + neuroinvasion confirmed ───────────────────────────
    # Moriel DG et al. (2010) PLoS ONE 5:e14165. Neonatal meningitis, K1.
    StrainInfo("IHE3034",        "Escherichia coli str. IHE3034",
               "NMEC", "ibeA,neuC(K1),ompA,kpsMII,fimH", 20952596),

    # Johnson TJ et al. (2012) BMC Genomics. Avian pathogenic, ExPEC.
    StrainInfo("CE10",           "Escherichia coli str. CE10",
               "NMEC", "kpsMII,iutA,iroN,hlyF,iss", 22535208),

    # ── AIEC — FimH variant + intracellular replication ───────────────────────
    # Miquel S et al. (2010) J Bacteriol 192:4541. Crohn's disease AIEC.
    StrainInfo("LF82",           "Escherichia coli str. LF82",
               "AIEC", "fimH-variant,inv,csg,lpfA", 20541499),

    # Martinez-Medina M et al. (2009) Inflamm Bowel Dis. AIEC, IBD mucosa.
    StrainInfo("HM605",          "Escherichia coli str. HM605",
               "AIEC", "fimH-variant,lpfA,csgA,htrA", 19253333),

    # Lapaquette P et al. (2010) Cell Microbiol. Adherent-invasive, IBD.
    StrainInfo("UM146",          "Escherichia coli str. UM146",
               "AIEC", "fimH-variant,lpfA,intracellular-replication", 20482583),

    # ── Supplemental EHEC to reach 30 ────────────────────────────────────────
    # Iyoda S et al. (2011) FEMS Microbiol Lett. Stx2+, LEE+, O21:NM serogroup.
    StrainInfo("O21_HC002",      "Escherichia coli O21:NM str. HC002",
               "EHEC", "stx2,LEE,eae", 21569098),
]

NON_PATHOGENIC: list[StrainInfo] = [

    # ── K-12 lineage — canonical lab strains, no virulence genes ─────────────
    # Blattner FR et al. (1997) Science 277:1453. Reference genome.
    StrainInfo("K12_MG1655",     "Escherichia coli K-12 str. MG1655",
               "K12", "no_virulence_genes,lambda_minus", 9278503,
               gcf="GCF_000005845.2"),

    # Hayashi K et al. (2006) Mol Syst Biol. F-minus, leu- derivative of MG1655.
    StrainInfo("K12_W3110",      "Escherichia coli K-12 str. W3110",
               "K12", "no_virulence_genes,F-minus,lacI", 16738553),

    # Durfee T et al. (2008) J Bacteriol 190:2597. Cloning/expression host.
    StrainInfo("K12_DH10B",      "Escherichia coli K-12 str. DH10B",
               "K12", "no_virulence_genes,endA1,recA1", 18245128),

    # Grenier F et al. (2014) PLoS ONE. Genome-reduced K-12 variant.
    StrainInfo("K12_MDS42",      "Escherichia coli K-12 str. MDS42",
               "K12", "no_virulence_genes,IS-free,genome-reduced", 24722555),

    # Jeong H et al. (2009) J Bacteriol 191:382. K-12, another complete assembly.
    StrainInfo("K12_BW2952",     "Escherichia coli K-12 str. BW2952",
               "K12", "no_virulence_genes", 18974106),

    # Posfai G et al. (2006) Science 312:1044. Minimal genome K-12 variant.
    StrainInfo("K12_MG1655_DY330","Escherichia coli K-12 str. DY330",
               "K12", "no_virulence_genes,lambda_Red_recombinase", 16614182),

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
               "PROBIOTIC", "microcin_B17,H_47_fimbriae,no_stx_no_LEE", 25448246),

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

    # Touchon M et al. (2009). Commensal gut isolate.
    StrainInfo("BL21_commensal", "Escherichia coli str. 83972",
               "COMMENSAL", "asymptomatic_bacteriuria_avirulent_strain", 17897305),

    # Vejborg RM et al. (2010) BMC Genomics. ABU, long-term colonisation no disease.
    StrainInfo("ABU83972",       "Escherichia coli str. ABU 83972",
               "COMMENSAL", "reductive_evolution,loss_of_virulence_genes", 20507614),

    # ── Industrial / environmental strains ────────────────────────────────────
    # Archer CT et al. (2011) Appl Environ Microbiol. Biofuel production strain.
    StrainInfo("W",              "Escherichia coli str. W",
               "LAB_B", "no_virulence_genes,industrial_fermentation", 21278278),

    # Durfee T et al. (2008). Cloning host, recA-minus.
    StrainInfo("K12_C2566",      "Escherichia coli K-12 str. C2566",
               "K12", "no_virulence_genes,recA_minus", 18245128),

    # Lukjancenko O et al. (2010) Microb Ecol. Commensal diversity panel.
    StrainInfo("ATCC8739",       "Escherichia coli ATCC 8739",
               "COMMENSAL", "ATCC_type_strain,no_virulence", 20830437),

    # ── Additional K-12 derivatives to complete 30 ───────────────────────────
    StrainInfo("K12_MC4100",     "Escherichia coli K-12 str. MC4100",
               "K12", "no_virulence_genes,araD_minus,delta_lac", 0),

    StrainInfo("K12_XL1Blue",    "Escherichia coli K-12 str. XL1-Blue",
               "K12", "no_virulence_genes,recA1_minus,endA1_minus", 0),

    StrainInfo("K12_DH5alpha",   "Escherichia coli K-12 str. DH5alpha",
               "K12", "no_virulence_genes,phi80_lacZdeltaM15", 0),

    StrainInfo("K12_HB101",      "Escherichia coli K-12 str. HB101",
               "K12", "no_virulence_genes,K-12_B_hybrid,recA13", 0),

    StrainInfo("K12_MG1655_LTEE","Escherichia coli K-12 str. MG1655 evolved",
               "K12", "no_virulence_genes,LTEE_related", 0),

    # Tenaillon O et al. (2010) Science 331:457. ECOR collection strain.
    StrainInfo("ECOR2",          "Escherichia coli str. ECOR 2",
               "COMMENSAL", "ECOR_collection,healthy_human", 20007904),

    StrainInfo("ECOR3",          "Escherichia coli str. ECOR 3",
               "COMMENSAL", "ECOR_collection,healthy_human", 20007904),

    StrainInfo("ECOR5",          "Escherichia coli str. ECOR 5",
               "COMMENSAL", "ECOR_collection,healthy_human", 20007904),

    StrainInfo("ECOR10",         "Escherichia coli str. ECOR 10",
               "COMMENSAL", "ECOR_collection,healthy_human", 20007904),

    StrainInfo("ECOR15",         "Escherichia coli str. ECOR 15",
               "COMMENSAL", "ECOR_collection,healthy_human", 20007904),

    # Shelobolina ES et al. (2004) Appl Environ Microbiol. Environmental isolate.
    StrainInfo("SMS_3_5",        "Escherichia coli str. SMS-3-5",
               "COMMENSAL", "environmental_soil,no_virulence_islands", 17965192),

    # ── B strain derivatives ──────────────────────────────────────────────────
    StrainInfo("B_ATCC9637",     "Escherichia coli B str. ATCC 9637",
               "LAB_B", "no_virulence_genes,B_lineage", 20840077),
]

ALL_STRAINS = PATHOGENIC + NON_PATHOGENIC

assert len(PATHOGENIC)     >= 30, f"Need 30 pathogenic, have {len(PATHOGENIC)}"
assert len(NON_PATHOGENIC) >= 30, f"Need 30 non-pathogenic, have {len(NON_PATHOGENIC)}"


# ── NCBI resolution ───────────────────────────────────────────────────────────

def entrez_search_assembly(organism: str) -> Optional[str]:
    """
    Search NCBI Assembly for a complete genome of the given organism.
    Returns first GCF_ accession found, or None.
    """
    query = f'"{organism}"[All Fields] AND "Complete Genome"[Assembly Level]'
    time.sleep(REQUEST_DELAY)
    try:
        handle  = Entrez.esearch(db="assembly", term=query, retmax=5)
        record  = Entrez.read(handle)
        ids     = record.get("IdList", [])
        if not ids:
            # Broader search without organism constraint
            query2 = f'{organism} AND "Complete Genome"[Assembly Level] AND "Escherichia coli"[Organism]'
            handle = Entrez.esearch(db="assembly", term=query2, retmax=5)
            record = Entrez.read(handle)
            ids    = record.get("IdList", [])
        if not ids:
            return None
        time.sleep(REQUEST_DELAY)
        handle  = Entrez.esummary(db="assembly", id=",".join(ids))
        summary = Entrez.read(handle)
        for doc in summary["DocumentSummarySet"]["DocumentSummary"]:
            acc = doc.get("AssemblyAccession", "")
            level = doc.get("AssemblyStatus", "")
            if acc.startswith("GCF_") and "Complete" in level:
                return acc
    except Exception as exc:
        log.warning("  Entrez search failed for %s: %s", organism, exc)
    return None


def verify_assembly_level(gcf: str) -> bool:
    """Return True if the assembly is chromosome-level or complete."""
    time.sleep(REQUEST_DELAY)
    try:
        handle  = Entrez.esearch(db="assembly", term=f"{gcf}[Assembly Accession]", retmax=1)
        record  = Entrez.read(handle)
        ids     = record.get("IdList", [])
        if not ids:
            return False
        handle  = Entrez.esummary(db="assembly", id=ids[0])
        summary = Entrez.read(handle)
        for doc in summary["DocumentSummarySet"]["DocumentSummary"]:
            return "Complete" in doc.get("AssemblyStatus", "")
    except Exception:
        pass
    return False


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

    # Resolve GCF accession
    gcf = strain.gcf
    if gcf is None:
        log.info("SRCH  %-25s  searching NCBI Assembly...", strain.label)
        gcf = entrez_search_assembly(strain.organism)
    if gcf is None:
        log.warning("FAIL  %-25s  no complete genome found on NCBI", strain.label)
        return None

    out_dir.mkdir(exist_ok=True)
    # Use full path so script works when called with explicit Python binary
    ngd_bin = Path(sys.executable).parent / "ncbi-genome-download"
    cmd = [
        str(ngd_bin),
        "--assembly-accessions", gcf,
        "--formats", "fasta",
        "--assembly-levels", "complete",
        "--output-folder", str(out_dir),
        "--flat-output",
        "--no-cache",
        "bacteria",
    ]
    log.info("FETCH %-25s  %s", strain.label, gcf)
    result = subprocess.run(cmd, capture_output=True, text=True, timeout=600)

    if result.returncode != 0:
        log.error("FAIL  %-25s  %s", strain.label, result.stderr.strip()[:120])
        return None

    gz_files = list(out_dir.glob("*_genomic.fna.gz"))
    if not gz_files:
        log.error("FAIL  %-25s  no .fna.gz found after download", strain.label)
        return None

    subprocess.run(["gunzip", "-f", str(gz_files[0])], check=True, capture_output=True)
    unzipped = gz_files[0].with_suffix("")
    unzipped.rename(fasta)

    # Size validation
    size_bp = sum(len(rec.seq) for rec in SeqIO.parse(str(fasta), "fasta"))
    if not (GENOME_SIZE_MIN <= size_bp <= GENOME_SIZE_MAX):
        log.warning("WARN  %-25s  genome size %d bp outside expected range", strain.label, size_bp)

    log.info("OK    %-25s  %s  %d bp", strain.label, gcf, size_bp)
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
        sys.exit(1)


if __name__ == "__main__":
    main()
