#!/usr/bin/env python3
"""
reblast_corrected_strains.py
Re-runs BLAST for the two corrected genomes (O103H2_12009, O26H11_11368)
against the existing 298 pilot markers. Prints corrected detection rates.
"""

import os
import subprocess
import tempfile
import pandas as pd

BLASTN     = "/opt/anaconda3/envs/my_ucd_env/bin/blastn"
MAKEBLAST  = "/opt/anaconda3/envs/my_ucd_env/bin/makeblastdb"

BASE       = "/Users/user/Documents/Projects/bioinformatics"
MARKERS    = f"{BASE}/data/results/pilot_heatmap/sakai_divergent_filtered.fna"
GENOME_DIR = f"{BASE}/data/genomes"

TARGETS = {
    "O103H2_12009": f"{GENOME_DIR}/O103H2_12009/O103H2_12009.fna",
    "O26H11_11368": f"{GENOME_DIR}/O26H11_11368/O26H11_11368.fna",
}

def count_markers(fna_path):
    n = 0
    with open(fna_path) as fh:
        for line in fh:
            if line.startswith(">"):
                n += 1
    return n

n_markers = count_markers(MARKERS)
print(f"Pilot marker set: {n_markers} markers")

with tempfile.TemporaryDirectory() as tmpdir:
    # Build combined FASTA with strain ID prefix
    combined = os.path.join(tmpdir, "targets.fna")
    with open(combined, "w") as out:
        for sid, fpath in TARGETS.items():
            print(f"  Adding {sid} ← {fpath}")
            with open(fpath) as fh:
                for line in fh:
                    if line.startswith(">"):
                        contig = line[1:].split()[0]
                        out.write(f">{sid}|{contig}\n")
                    else:
                        out.write(line)

    # Build BLAST database
    db = os.path.join(tmpdir, "targets_db")
    r = subprocess.run(
        [MAKEBLAST, "-in", combined, "-dbtype", "nucl", "-out", db, "-parse_seqids"],
        capture_output=True, text=True
    )
    if r.returncode != 0:
        raise RuntimeError(f"makeblastdb failed:\n{r.stderr[:400]}")
    print("  BLAST database built")

    # Run BLASTn
    blast_out = os.path.join(tmpdir, "blast.tsv")
    cols = "qseqid sseqid pident length qlen qstart qend sstart send evalue bitscore"
    r = subprocess.run(
        [BLASTN, "-query", MARKERS, "-db", db,
         "-outfmt", f"6 {cols}",
         "-perc_identity", "85", "-qcov_hsp_perc", "50",
         "-max_target_seqs", "500", "-num_threads", "4",
         "-out", blast_out],
        capture_output=True, text=True
    )
    if r.returncode != 0:
        raise RuntimeError(f"blastn failed:\n{r.stderr[:400]}")

    # Parse results
    with open(blast_out) as fh:
        nhits = sum(1 for _ in fh)
    print(f"  Total BLAST hits: {nhits}")

    df = pd.read_csv(blast_out, sep="\t", header=None, names=cols.split())
    df["strain"] = df["sseqid"].str.split("|").str[0]
    df["qcov"] = (df["qend"] - df["qstart"] + 1) / df["qlen"]
    df = df[df["qcov"] >= 0.5]

    print(f"\n{'Strain':<22} {'Markers present':>15}  {'Pct':>6}")
    print("-" * 46)
    for sid in TARGETS:
        present = df[df["strain"] == sid]["qseqid"].nunique()
        pct = present / n_markers * 100
        print(f"{sid:<22} {present:>8}/{n_markers:<6}  {pct:>5.1f}%")

print("\nDone.")
