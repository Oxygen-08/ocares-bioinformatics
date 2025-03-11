#!/usr/bin/env python3
"""
Extract Divergent Regions from BLAST Results

This script:
1. Reads BLAST results and identifies unique divergent regions
2. Extracts sequences from the original genome files
3. Prepares sequences for functional annotation
4. Handles both nucleotide and protein sequence extraction (with translation option)
5. Includes proper sequence headers with coordinates and genome information

Usage:
    python extract_divergent_regions.py \
        --blast_results results/blast_results.txt \
        --genomes_dir genomes/ \
        --output_dir results/divergent_regions/ \
        --min_identity 90 \
        --max_identity 100 \
        --min_length 50 \
        --translate False \
        --genome_mapping results/adapted_genome_name_mapping.csv
"""

import os
import argparse
import pandas as pd
import numpy as np
from collections import defaultdict
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import re
import logging

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger('extract_divergent_regions')

def parse_args():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description='Extract divergent regions from BLAST results')
    
    parser.add_argument('--blast_results', type=str, required=True,
                        help='Path to BLAST results file in tabular format')
    parser.add_argument('--genomes_dir', type=str, required=True,
                        help='Directory containing genome FASTA files')
    parser.add_argument('--output_dir', type=str, default='results/divergent_regions',
                        help='Output directory for extracted sequences')
    parser.add_argument('--min_identity', type=float, default=0.0,
                        help='Minimum percent identity for regions to include')
    parser.add_argument('--max_identity', type=float, default=95.0,
                        help='Maximum percent identity for regions to be considered divergent')
    parser.add_argument('--min_length', type=int, default=50,
                        help='Minimum alignment length to consider')
    parser.add_argument('--translate', type=bool, default=False,
                        help='Translate nucleotide sequences to protein')
    parser.add_argument('--genome_mapping', type=str, default=None,
                        help='Path to genome name mapping file')
    
    return parser.parse_args()

def read_blast_results(blast_file):
    """Read BLAST results from tabular format."""
    try:
        column_names = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 
                        'gapopen', 'qstart', 'qend', 'sstart', 'send', 
                        'evalue', 'bitscore']
        
        df = pd.read_csv(blast_file, sep='\t', header=None, names=column_names)
        
        logger.info(f"Successfully read {len(df)} BLAST alignments from {blast_file}")
        return df
    except Exception as e:
        logger.error(f"Error reading BLAST results: {e}")
        raise

def load_genome_mapping(mapping_file):
    """Load genome name mapping if provided."""
    if mapping_file and os.path.exists(mapping_file):
        try:
            mapping_df = pd.read_csv(mapping_file)
            # Create dictionary with Original_ID as key and Simple_Name as value
            mapping_dict = dict(zip(mapping_df['Original_ID'], mapping_df['Simple_Name']))
            logger.info(f"Loaded genome mapping with {len(mapping_dict)} entries")
            return mapping_dict
        except Exception as e:
            logger.warning(f"Failed to load genome mapping: {e}")
    
    return {}

def extract_region_info(region_id):
    """Extract genome and coordinate information from region ID."""
    # Match pattern like NZ_CP101522.1:12345-67890
    match = re.match(r'([^:]+):(\d+)-(\d+)', region_id)
    
    if match:
        genome_id = match.group(1)
        start = int(match.group(2))
        end = int(match.group(3))
        return genome_id, start, end
    
    # Handle other formats if present
    logger.warning(f"Could not parse region ID format: {region_id}")
    return region_id, None, None

def find_genome_file(genome_id, genomes_dir):
    """Find the genome file corresponding to the genome ID."""
    # Try exact match first
    for file in os.listdir(genomes_dir):
        if file.endswith(('.fa', '.fasta', '.fna')) and genome_id in file:
            return os.path.join(genomes_dir, file)
    
    # Try a more flexible match if exact match fails
    base_id = genome_id.split('.')[0]
    for file in os.listdir(genomes_dir):
        if file.endswith(('.fa', '.fasta', '.fna')) and base_id in file:
            return os.path.join(genomes_dir, file)
    
    logger.warning(f"Could not find genome file for ID: {genome_id}")
    return None

def extract_sequence(genome_file, start, end, reverse_complement=False):
    """Extract a sequence from a genome file using coordinates."""
    try:
        for record in SeqIO.parse(genome_file, "fasta"):
            # Get the sequence
            seq = record.seq[start-1:end]
            
            if reverse_complement:
                seq = seq.reverse_complement()
                
            return seq
    except Exception as e:
        logger.error(f"Error extracting sequence: {e}")
        return None

def translate_sequence(seq):
    """Translate nucleotide sequence to protein."""
    try:
        # Find the best ORF by translating in all 6 reading frames
        translations = []
        
        # Forward frames
        for i in range(3):
            trans = seq[i:].translate(to_stop=True)
            if len(trans) > 0:
                translations.append((len(trans), trans, i, 1))
        
        # Reverse frames
        rev_seq = seq.reverse_complement()
        for i in range(3):
            trans = rev_seq[i:].translate(to_stop=True)
            if len(trans) > 0:
                translations.append((len(trans), trans, i, -1))
        
        # Sort by length, get the longest
        if translations:
            translations.sort(reverse=True)
            return translations[0][1], translations[0][2], translations[0][3]
        
        # If no ORF found, return the translated sequence in frame 0
        return seq.translate(), 0, 1
    except Exception as e:
        logger.error(f"Error translating sequence: {e}")
        return None, None, None

def identify_divergent_regions(blast_df, min_identity, max_identity, min_length):
    """Identify divergent regions based on percent identity and length filters."""
    # Filter alignments based on criteria
    filtered_df = blast_df[
        (blast_df['pident'] >= min_identity) & 
        (blast_df['pident'] <= max_identity) & 
        (blast_df['length'] >= min_length)
    ]
    
    # Group by query identifier
    divergent_regions = defaultdict(list)
    for _, row in filtered_df.iterrows():
        region_id = row['qseqid']
        divergent_regions[region_id].append(row)
    
    logger.info(f"Identified {len(divergent_regions)} divergent regions")
    return divergent_regions

def extract_and_save_sequences(divergent_regions, genomes_dir, output_dir, 
                               translate=False, genome_mapping=None):
    """Extract and save sequences for divergent regions."""
    os.makedirs(output_dir, exist_ok=True)
    
    # Create output files
    nucleotide_file = os.path.join(output_dir, "divergent_regions_nucleotide.fasta")
    protein_file = os.path.join(output_dir, "divergent_regions_protein.fasta") if translate else None
    
    # Open output files
    with open(nucleotide_file, 'w') as nuc_out:
        prot_out = open(protein_file, 'w') if translate else None
        
        extracted_count = 0
        failed_count = 0
        
        # Process each divergent region
        for region_id, alignments in divergent_regions.items():
            genome_id, start, end = extract_region_info(region_id)
            
            if not all([genome_id, start, end]):
                logger.warning(f"Skipping region with incomplete coordinates: {region_id}")
                failed_count += 1
                continue
            
            # Get genome filename
            genome_file = find_genome_file(genome_id, genomes_dir)
            if not genome_file:
                failed_count += 1
                continue
            
            # Extract sequence
            sequence = extract_sequence(genome_file, start, end)
            if not sequence:
                failed_count += 1
                continue
            
            # Format display name using mapping if available
            display_name = genome_mapping.get(region_id, genome_id) if genome_mapping else genome_id
            
            # Create header with detailed information
            header = f"{display_name}|{start}-{end}|len={end-start+1}"
            
            # Write nucleotide sequence
            nuc_out.write(f">{header}\n{sequence}\n")
            
            # Translate if requested
            if translate and prot_out:
                protein, frame, strand = translate_sequence(sequence)
                if protein:
                    strand_symbol = "+" if strand > 0 else "-"
                    prot_header = f"{header}|frame={frame}|strand={strand_symbol}"
                    prot_out.write(f">{prot_header}\n{protein}\n")
            
            extracted_count += 1
        
        if prot_out:
            prot_out.close()
    
    logger.info(f"Extracted {extracted_count} sequences, failed to extract {failed_count}")
    logger.info(f"Nucleotide sequences saved to: {nucleotide_file}")
    if translate:
        logger.info(f"Protein sequences saved to: {protein_file}")
    
    return nucleotide_file, protein_file

def prepare_annotation_inputs(nucleotide_file, protein_file, output_dir):
    """Prepare files for functional annotation software."""
    annotation_dir = os.path.join(output_dir, "annotation_input")
    os.makedirs(annotation_dir, exist_ok=True)
    
    # Files for different annotation tools
    prokka_file = os.path.join(annotation_dir, "for_prokka.fasta")
    interproscan_file = os.path.join(annotation_dir, "for_interproscan.fasta")
    
    # For tools that use nucleotide input
    if nucleotide_file and os.path.exists(nucleotide_file):
        nuc_records = list(SeqIO.parse(nucleotide_file, "fasta"))
        
        # For Prokka
        with open(prokka_file, 'w') as f_out:
            SeqIO.write(nuc_records, f_out, "fasta")
        logger.info(f"Prepared nucleotide sequences for Prokka: {prokka_file}")
    
    # For tools that use protein input
    if protein_file and os.path.exists(protein_file):
        prot_records = list(SeqIO.parse(protein_file, "fasta"))
        
        # For InterProScan
        with open(interproscan_file, 'w') as f_out:
            SeqIO.write(prot_records, f_out, "fasta")
        logger.info(f"Prepared protein sequences for InterProScan: {interproscan_file}")
    
    # Return paths to annotation input files
    return {
        'prokka': prokka_file if os.path.exists(prokka_file) else None,
        'interproscan': interproscan_file if os.path.exists(interproscan_file) else None
    }

def create_annotation_commands(annotation_inputs, output_dir):
    """Create command templates for common annotation tools."""
    command_file = os.path.join(output_dir, "annotation_commands.sh")
    
    with open(command_file, 'w') as f:
        f.write("#!/bin/bash\n\n")
        f.write("# Commands for functional annotation of divergent regions\n\n")
        
        # Prokka command
        if annotation_inputs.get('prokka'):
            f.write("# Prokka annotation (bacterial gene prediction and annotation)\n")
            f.write(f"prokka --outdir {output_dir}/prokka_results \\\n")
            f.write(f"       --prefix divergent_regions \\\n")
            f.write(f"       --locustag DIV \\\n")
            f.write(f"       --compliant \\\n")
            f.write(f"       {annotation_inputs['prokka']}\n\n")
        
        # InterProScan command
        if annotation_inputs.get('interproscan'):
            f.write("# InterProScan annotation (protein domain and family prediction)\n")
            f.write(f"interproscan.sh -i {annotation_inputs['interproscan']} \\\n")
            f.write(f"               -f TSV,GFF3,HTML \\\n")
            f.write(f"               -d {output_dir}/interproscan_results \\\n")
            f.write(f"               -goterms -iprlookup -pa\n\n")
        
        # BLAST against nr/nt command
        if annotation_inputs.get('prokka'):
            f.write("# BLAST against NCBI nr/nt database\n")
            f.write(f"blastn -query {annotation_inputs['prokka']} \\\n")
            f.write(f"       -db nt \\\n")
            f.write(f"       -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle' \\\n")
            f.write(f"       -evalue 1e-5 \\\n")
            f.write(f"       -max_target_seqs 5 \\\n")
            f.write(f"       -out {output_dir}/blast_nr_results.txt\n\n")
        
        # HMMER command for protein family annotation
        if annotation_inputs.get('interproscan'):
            f.write("# HMMER search against Pfam database\n")
            f.write(f"hmmscan --tblout {output_dir}/hmmer_pfam_results.txt \\\n")
            f.write(f"        --noali \\\n")
            f.write

