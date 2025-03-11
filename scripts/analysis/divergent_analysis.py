import subprocess
import os

def identify_divergent_regions(reference, query, output_prefix):
    """
    Function to identify divergent regions in the reference genome
    by aligning against the query genome using MUMmer's nucmer.
    
    Args:
    - reference (str): Path to the reference genome FASTA file.
    - query (str): Path to the query genome FASTA file.
    - output_prefix (str): Prefix for output files.
    
    Outputs:
    - <output_prefix>.coords: Alignment coordinates file.
    - aligned_regions.txt: Aligned regions in reference.
    - divergent_regions.bed: Regions unique to the reference.
    """
    
    try:
        # Step 1: Run nucmer with --mum to find non-overlapping matches
        nucmer_cmd = f"nucmer --mum {reference} {query} --prefix={output_prefix}"
        subprocess.run(nucmer_cmd, shell=True, check=True)
        print("Step 1: NUCmer alignment completed.")

        # Step 2: Extract alignment coordinates
        coords_file = f"{output_prefix}.coords"
        show_coords_cmd = f"show-coords -rcl {output_prefix}.delta > {coords_file}"
        subprocess.run(show_coords_cmd, shell=True, check=True)
        print("Step 2: Extracted alignment coordinates.")

        # Step 3: Extract reference-aligned regions (start, end positions)
        aligned_regions_file = "aligned_regions.txt"
        awk_cmd = f"awk '{{print $1, $2}}' {coords_file} > {aligned_regions_file}"
        subprocess.run(awk_cmd, shell=True, check=True)
        print("Step 3: Extracted aligned regions.")

        # Step 4: Index reference genome for bedtools complement
        samtools_cmd = f"samtools faidx {reference}"
        subprocess.run(samtools_cmd, shell=True, check=True)
        print("Step 4: Indexed reference genome.")

        # Step 5: Identify divergent regions (regions not aligned)
        fai_file = f"{reference}.fai"
        divergent_regions_file = "divergent_regions.bed"
        bedtools_cmd = f"bedtools complement -i {aligned_regions_file} -g {fai_file} > {divergent_regions_file}"
        subprocess.run(bedtools_cmd, shell=True, check=True)
        print("Step 5: Identified divergent regions.")

        print("Analysis completed successfully!")

    except subprocess.CalledProcessError as e:
        print(f"Error during execution: {e}")
