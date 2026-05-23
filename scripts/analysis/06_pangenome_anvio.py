#!/usr/bin/env python3
"""
Anvi'o pangenome construction — must run inside the anvio8 conda environment.

    conda activate anvio8
    python scripts/analysis/06_pangenome_anvio.py

Or without activating:
    conda run -n anvio8 python scripts/analysis/06_pangenome_anvio.py

The genome manifest (data/genomes/genome_manifest.tsv) is the single source
of truth for which genomes are included. Any reformatted FASTA or contigs
database that does not correspond to a manifest entry is removed.

PATHOGEN/COMMENSAL classification:
    PATHOGEN  — EHEC, EPEC, UPEC, ETEC, NMEC, EAEC, AIEC, APEC
    COMMENSAL — K12, LAB_B, COMMENSAL, PROBIOTIC

Use --force to delete existing GENOMES.db and PAN.db and rerun from scratch.
"""

import argparse
import csv
import datetime
import os
import shutil
import subprocess
import sys
from pathlib import Path

ANVIO_BIN = Path("/opt/anaconda3/envs/anvio8/bin")
THREADS   = "4"
PROJECT   = "EcoliPangenome"

PATHOGENIC_PATHOTYPES = {"EHEC", "EPEC", "UPEC", "ETEC", "NMEC", "EAEC", "AIEC", "APEC"}

REPO_ROOT     = Path(__file__).resolve().parents[2]
GENOME_DIR    = REPO_ROOT / "data" / "genomes"
PAN_DIR       = REPO_ROOT / "data" / "results" / "pangenome"
REFORM_DIR    = PAN_DIR / "reformatted"
DB_DIR        = PAN_DIR / "contigs_dbs"
MANIFEST      = GENOME_DIR / "genome_manifest.tsv"
EXT_GENOMES   = PAN_DIR / "external_genomes.txt"
GENOME_LIST   = PAN_DIR / "genome_list.txt"
LAYER_GROUPS  = PAN_DIR / "layer_groups.txt"
ENRICHMENT    = PAN_DIR / "enrichment_scores.tsv"
GENE_CLUSTERS = PAN_DIR / "gene_clusters.tsv"
STORAGE       = PAN_DIR / f"{PROJECT}-GENOMES.db"
PAN_DB        = PAN_DIR / f"{PROJECT}-PAN.db"


def ts():
    return datetime.datetime.now().strftime("%H:%M:%S")


def log(msg):
    print(f"[{ts()}] {msg}", flush=True)


def run(cmd, desc="", check=True):
    env = {**os.environ, "PATH": f"{ANVIO_BIN}:{os.environ['PATH']}"}
    label = desc or Path(cmd[0]).name
    log(f"  → {label}")
    result = subprocess.run([str(c) for c in cmd], env=env,
                            capture_output=True, text=True)
    if check and result.returncode != 0:
        print(result.stderr[-800:], file=sys.stderr)
        raise RuntimeError(f"Failed: {label}")
    return result


def safe_label(label: str) -> str:
    """Anvi'o names must start with a letter."""
    return f"G_{label}" if label[0].isdigit() else label


def orig_label(slabel: str) -> str:
    return slabel[2:] if slabel.startswith("G_") else slabel


# ── Parse args ────────────────────────────────────────────────────────────────
parser = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("--force", action="store_true",
                    help="Delete existing GENOMES.db and PAN.db and rebuild from scratch")
args = parser.parse_args()

PAN_DIR.mkdir(parents=True, exist_ok=True)
REFORM_DIR.mkdir(parents=True, exist_ok=True)
DB_DIR.mkdir(parents=True, exist_ok=True)


# ── 1. Read manifest ──────────────────────────────────────────────────────────
log("Reading genome manifest")
manifest = {}
with open(MANIFEST) as f:
    for row in csv.DictReader(f, delimiter="\t"):
        manifest[row["label"]] = row["pathotype"]

pathogen_set  = {l for l, p in manifest.items() if p in PATHOGENIC_PATHOTYPES}
commensal_set = {l for l, p in manifest.items() if p not in PATHOGENIC_PATHOTYPES}
log(f"  {len(manifest)} genomes | PATHOGEN={len(pathogen_set)} COMMENSAL={len(commensal_set)}")


# ── 2. Remove stale files not in manifest ─────────────────────────────────────
log("Purging stale files not in manifest")
for fna in list(REFORM_DIR.glob("*.fna")):
    if fna.stem not in manifest:
        log(f"  remove stale reformatted: {fna.name}")
        fna.unlink()
for db in list(DB_DIR.glob("*.db")):
    if db.stem not in manifest:
        log(f"  remove stale contigs_db: {db.name}")
        db.unlink()


# ── 3. Force-delete upstream DBs and BLAST intermediates if requested ──────────
if args.force:
    log("--force: removing GENOMES.db, PAN.db, BLAST intermediates, ANI")
    blast_intermediates = [
        "blast-search-results.txt",
        "combined-aas.fa",
        "combined-aas.fa.unique",
        "combined-aas.fa.unique.names",
        "combined-aas.fa.unique.pdb",
        "combined-aas.fa.unique.phr",
        "combined-aas.fa.unique.pin",
        "combined-aas.fa.unique.pjs",
        "combined-aas.fa.unique.pot",
        "combined-aas.fa.unique.psq",
        "combined-aas.fa.unique.ptf",
        "combined-aas.fa.unique.pto",
        "mcl-input.txt",
        "mcl-clusters.txt",
        "log.txt",
    ]
    for name in blast_intermediates:
        p = PAN_DIR / name
        if p.exists():
            p.unlink()
            log(f"  removed {name}")
    for p in [STORAGE, PAN_DB, PAN_DIR / "ANI"]:
        if p.exists():
            shutil.rmtree(p) if p.is_dir() else p.unlink()
            log(f"  removed {p.name}")


# ── 4. Reformat + build contigs_db + annotate for each manifest genome ─────────
log("Processing manifest genomes")
for label, pathotype in sorted(manifest.items()):
    fna_src = GENOME_DIR / label / f"{label}.fna"
    fna_out = REFORM_DIR / f"{label}.fna"
    db_path = DB_DIR / f"{label}.db"
    slabel  = safe_label(label)

    if not fna_src.exists():
        log(f"  WARNING: {label}.fna not found — skipping")
        continue

    if not fna_out.exists():
        run([ANVIO_BIN / "anvi-script-reformat-fasta",
             fna_src, "-o", fna_out,
             "--min-len", "1000",
             "--simplify-names",
             "--prefix", slabel,
             "--seq-type", "NT"],
            desc=f"reformat {label}")

    if not db_path.exists():
        run([ANVIO_BIN / "anvi-gen-contigs-database",
             "-f", fna_out, "-o", db_path,
             "-n", slabel, "--num-threads", THREADS],
            desc=f"contigs_db {label}")

        run([ANVIO_BIN / "anvi-run-hmms",
             "-c", db_path,
             "-I", "Bacteria_71",
             "--num-threads", THREADS, "--quiet"],
            desc=f"HMMs {label}")

        run([ANVIO_BIN / "anvi-run-ncbi-cogs",
             "-c", db_path,
             "--cog-version", "COG14",
             "--num-threads", THREADS],
            desc=f"COG14 {label}")


# ── 5. Write external_genomes.txt and genome_list.txt ─────────────────────────
log("Writing external_genomes.txt and genome_list.txt")
entries = []
for label in sorted(manifest.keys()):
    db_path = DB_DIR / f"{label}.db"
    if db_path.exists():
        entries.append((safe_label(label), label, db_path))

with open(EXT_GENOMES, "w") as f:
    f.write("name\tcontigs_db_path\n")
    for slabel, _, db_path in entries:
        f.write(f"{slabel}\t{db_path}\n")

with open(GENOME_LIST, "w") as f:
    for slabel, label, _ in entries:
        f.write(f"{slabel}\t{REFORM_DIR / (label + '.fna')}\n")

log(f"  {len(entries)} genomes registered")


# ── 6. Write layer_groups.txt ─────────────────────────────────────────────────
log("Writing layer_groups.txt")
with open(LAYER_GROUPS, "w") as f:
    f.write("name\tgroup\n")
    for slabel, label, db_path in entries:
        group = "PATHOGEN" if label in pathogen_set else "COMMENSAL"
        f.write(f"{slabel}\t{group}\n")


# ── 7. Build genomes storage ──────────────────────────────────────────────────
if not STORAGE.exists():
    log("Building genomes storage")
    run([ANVIO_BIN / "anvi-gen-genomes-storage",
         "-e", EXT_GENOMES, "-o", STORAGE],
        desc="anvi-gen-genomes-storage")
else:
    log("Genomes storage exists — skipping (use --force to rebuild)")


# ── 8. Pangenome analysis ──────────────────────────────────────────────────────
if not PAN_DB.exists():
    log("Running pangenome (MCL inflation=10, minbit=0.8) — this takes 2–6 hours")
    run([ANVIO_BIN / "anvi-pan-genome",
         "-g", STORAGE,
         "--project-name", PROJECT,
         "--output-dir", PAN_DIR,
         "--num-threads", THREADS,
         "--minbit", "0.8",
         "--mcl-inflation", "10",
         "--use-ncbi-blast"],
        desc="anvi-pan-genome")
else:
    log("PAN.db exists — skipping (use --force to rebuild)")


# ── 9. Import layer groups into pan db ────────────────────────────────────────
log("Importing layer groups into pan db")
run([ANVIO_BIN / "anvi-import-misc-data",
     "-p", PAN_DB,
     "-t", "layers",
     LAYER_GROUPS],
    desc="import layer_groups",
    check=False)


# ── 10. Functional enrichment (COG14_FUNCTION) ────────────────────────────────
log("Computing functional enrichment (COG14_FUNCTION)")
run([ANVIO_BIN / "anvi-compute-functional-enrichment-in-pan",
     "-p", PAN_DB,
     "-g", STORAGE,
     "--category-variable", "group",
     "--annotation-source", "COG14_FUNCTION",
     "-o", ENRICHMENT],
    desc="anvi-compute-functional-enrichment-in-pan")
log(f"  enrichment_scores.tsv written")


# ── 11. Export gene cluster presence/absence matrix ───────────────────────────
log("Exporting gene_clusters.tsv")
run([ANVIO_BIN / "anvi-export-table",
     PAN_DB,
     "--table", "gene_clusters",
     "-o", GENE_CLUSTERS],
    desc="anvi-export-table gene_clusters")
log(f"  gene_clusters.tsv written")

log("=== Pangenome complete. Run 07_ml_classifier.py next. ===")
