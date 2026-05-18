"""
Rebuild data/metagenome/read_labels.tsv with correct is_pathogen labels.

InSilicoSeq uses NCBI accession numbers from FASTA headers as read ID prefixes
(e.g., NC_002695.2_0_0/1). The original build_read_label_table() in
04_simulate_metagenome.py matched against genome labels (e.g., O157H7_Sakai)
which never appear in FASTQ headers — so all 750K reads got is_pathogen=0.

Fix: build accession→(genome_label, pathotype) lookup from FASTA headers,
then label each read by accession prefix.
"""

import gzip
import re
from pathlib import Path

PROJ = Path("/Users/user/Documents/Projects/bioinformatics")
MANIFEST = PROJ / "data/genomes/genome_manifest.tsv"
META_DIR = PROJ / "data/metagenome"
OUT = META_DIR / "read_labels.tsv"

EHEC_PATHOTYPES = {"EHEC"}  # is_pathogen=1 for these

CONDITIONS = {
    "low":  META_DIR / "low/low_merged.fastq",
    "mid":  META_DIR / "mid/mid_merged.fastq",
    "high": META_DIR / "high/high_merged.fastq",
}


def extract_accession(read_id: str) -> str:
    """
    Convert read ID like 'NC_002695.2_0_0/1' → accession 'NC_002695.2'.
    InSilicoSeq format: <contig_name>_<position>_<frag>/[12]
    Contig names are 'NC_XXXXXX.X' or 'NZ_XXXXXX.X' — always two underscore-split
    fields for the accession (NC + 002695.2 or NZ + CP171242.1).
    """
    # Strip trailing /1 or /2
    rid = read_id.rstrip("/12").rstrip("/")
    # Handle /1 /2 suffix properly
    rid = re.sub(r'/[12]$', '', read_id)
    parts = rid.split("_")
    if len(parts) >= 2 and parts[0] in ("NC", "NZ"):
        return f"{parts[0]}_{parts[1]}"
    # Fallback: contig with no underscore (e.g., CP001164.1_0_0)
    return parts[0] if parts else read_id


def build_accession_lookup() -> dict:
    """Read genome manifest, parse each FASTA, return {accession: (label, pathotype)}."""
    lookup = {}
    with open(MANIFEST) as f:
        header = f.readline().strip().split("\t")
        label_idx = header.index("label")
        path_idx = header.index("pathotype")
        fasta_idx = header.index("fasta_path")
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) <= fasta_idx:
                continue
            genome_label = parts[label_idx]
            pathotype = parts[path_idx]
            fasta_path = Path(parts[fasta_idx])
            if not fasta_path.exists():
                print(f"  [WARN] FASTA not found: {fasta_path}")
                continue
            opener = gzip.open if fasta_path.suffix == ".gz" else open
            with opener(fasta_path, "rt") as fa:
                for line in fa:
                    if line.startswith(">"):
                        # ">NC_002695.2 Escherichia coli..." → "NC_002695.2"
                        contig_name = line[1:].split()[0]
                        accession = extract_accession(contig_name + "_0_0")
                        lookup[accession] = (genome_label, pathotype)
    print(f"  Built lookup: {len(lookup)} contig accessions from {MANIFEST.name}")
    return lookup


def rebuild_labels(lookup: dict) -> None:
    rows = []
    for condition, fastq_path in sorted(CONDITIONS.items()):
        if not fastq_path.exists():
            print(f"  [WARN] FASTQ not found: {fastq_path}")
            continue
        print(f"  Processing {condition} ({fastq_path.name})...")
        n_reads = 0
        n_matched = 0
        with open(fastq_path) as fq:
            while True:
                header_line = fq.readline()
                if not header_line:
                    break
                seq_line   = fq.readline()
                plus_line  = fq.readline()
                qual_line  = fq.readline()
                if not header_line.startswith("@"):
                    continue
                read_id = header_line[1:].strip()
                accession = extract_accession(read_id)
                if accession in lookup:
                    genome_label, pathotype = lookup[accession]
                    is_pathogen = 1 if pathotype in EHEC_PATHOTYPES else 0
                    n_matched += 1
                else:
                    genome_label = "unknown"
                    is_pathogen = 0
                rows.append(f"{read_id}\t{condition}\t{genome_label}\t{is_pathogen}\n")
                n_reads += 1
        pct = 100 * n_matched / n_reads if n_reads else 0
        print(f"    {n_reads:,} reads | {n_matched:,} matched ({pct:.1f}%)")

    with open(OUT, "w") as out:
        out.write("read_id\tcondition\tgenome\tis_pathogen\n")
        out.writelines(rows)
    n_pathogen = sum(1 for r in rows if r.endswith("\t1\n"))
    print(f"\n  Written: {OUT}")
    print(f"  Total reads: {len(rows):,} | Pathogen reads: {n_pathogen:,} ({100*n_pathogen/len(rows):.2f}%)")


if __name__ == "__main__":
    print("Building accession → genome lookup...")
    lookup = build_accession_lookup()
    print("\nLabelling reads...")
    rebuild_labels(lookup)
    print("\nDone.")
