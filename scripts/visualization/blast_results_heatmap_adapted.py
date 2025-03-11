#!/usr/bin/env python3
"""
Generate heatmap visualizations from BLAST results data specifically for genomes with accession
number format NZ_CP101522.1:XXXXX-XXXXX.

This adapted script:
1. Extracts coordinates from query regions in NZ_CP101522.1:XXXXX-XXXXX format
2. Extracts just the accession numbers from subject genomes without database prefixes
3. Creates optimized heatmaps with rotated labels for better readability

Original script modified to handle specific identifier formats used in the current project.
"""

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import scipy.cluster.hierarchy as hierarchy
from scipy.spatial.distance import pdist
import matplotlib.gridspec as gridspec
import argparse

def parse_blast_results(file_path):
    """
    Parse the BLAST results file and extract relevant information.
    
    Args:
        file_path (str): Path to the BLAST results file
        
    Returns:
        pandas.DataFrame: DataFrame containing parsed BLAST data
    """
    # Read in the BLAST results file
    # Assuming tab-separated format with standard BLAST output columns
    try:
        # Try with standard BLAST outfmt 6 column names
        columns = ['query_id', 'subject_id', 'percent_identity', 'alignment_length',
                  'mismatches', 'gap_opens', 'q_start', 'q_end', 's_start', 's_end',
                  'e_value', 'bit_score']
        
        df = pd.read_csv(file_path, sep='\t', header=None)
        
        # If the file has a header row, read it again with header
        if not df.iloc[0, 0].replace('.', '', 1).isdigit():
            df = pd.read_csv(file_path, sep='\t')
        else:
            # If no header, assign column names
            if df.shape[1] == len(columns):
                df.columns = columns
            else:
                # If column count doesn't match, use generic column names
                df.columns = [f'col_{i}' for i in range(df.shape[1])]
                
                # Try to identify key columns by content
                for i, col in enumerate(df.columns):
                    # Percent identity is usually a value between 0-100
                    if df[col].dtype in [np.float64, np.int64] and df[col].mean() > 0 and df[col].mean() <= 100:
                        df.rename(columns={col: 'percent_identity'}, inplace=True)
                        break
        
        # Extract query and subject identifiers if available
        if 'query_id' not in df.columns and df.shape[1] > 0:
            df.rename(columns={df.columns[0]: 'query_id'}, inplace=True)
        if 'subject_id' not in df.columns and df.shape[1] > 1:
            df.rename(columns={df.columns[1]: 'subject_id'}, inplace=True)
        if 'percent_identity' not in df.columns and df.shape[1] > 2:
            df.rename(columns={df.columns[2]: 'percent_identity'}, inplace=True)
            
        return df
    
    except Exception as e:
        print(f"Error parsing BLAST results: {e}")
        return None

def extract_region_info(query_id):
    """
    Extract coordinates from query ID in NZ_CP101522.1:XXXXX-XXXXX format.
    
    Args:
        query_id (str): Query ID from BLAST results
        
    Returns:
        str: Extracted coordinates (XXXXX-XXXXX)
    """
    query_id = str(query_id)
    
    # Handle format like NZ_CP101522.1:12345-67890
    if ':' in query_id:
        # Extract the part after the colon (coordinates)
        return query_id.split(':', 1)[1]
    
    # Handle format with pipe delimiter
    parts = query_id.split('|')
    if len(parts) > 1:
        for part in parts:
            if '-' in part:  # Look for coordinate pattern
                return part
    
    # If no specific format is detected, return the original ID
    return query_id

def create_genome_name_mapping(blast_df):
    """
    Create a mapping from original genome identifiers to simplified names.
    
    Args:
        blast_df (pandas.DataFrame): DataFrame containing parsed BLAST data
        
    Returns:
        dict: Mapping from original genome IDs to simplified names
    """
    # Extract all unique genome identifiers
    original_ids = []
    for subject_id in blast_df['subject_id'].unique():
        # Use the extract_genome_info function to get the cleaned identifier
        cleaned_id = extract_genome_info_original(subject_id)
        original_ids.append((subject_id, cleaned_id))
    
    # Create a mapping from original IDs to simple names
    mapping = {}
    for i, (original_id, cleaned_id) in enumerate(sorted(original_ids, key=lambda x: x[1])):
        simple_name = f"Genome_{i+1}"
        mapping[original_id] = simple_name
    
    return mapping

def extract_genome_info_original(subject_id):
    """
    Extract accession number from subject ID without database prefixes.
    Original version that doesn't use the simplified names.
    
    Args:
        subject_id (str): Subject ID from BLAST results
        
    Returns:
        str: Extracted accession number without database prefixes
    """
    subject_id = str(subject_id)
    
    # Handle NCBI gi format: gi|12345|ref|NC_123456.1|
    if subject_id.startswith('gi|'):
        # Find the accession number after 'ref|' if it exists
        if '|ref|' in subject_id:
            accession = subject_id.split('|ref|')[1].split('|')[0]
            return accession
    
    # Handle NZ_, NC_ format
    for prefix in ['NZ_', 'NC_']:
        if prefix in subject_id:
            # Extract just the accession part
            parts = subject_id.split(prefix, 1)
            if len(parts) > 1:
                # Return with the prefix (e.g., CP101522.1 from NZ_CP101522.1)
                accession = parts[1].split()[0].split('|')[0].split(':')[0]
                return accession
    
    # Handle formats with pipe delimiter
    parts = subject_id.split('|')
    if len(parts) > 1:
        # Try to find the part that looks like an accession number
        for part in parts:
            if any(part.startswith(prefix) for prefix in ['NZ_', 'NC_', 'CP', 'NW_']):
                return part.split(':')[0]
        # If no recognizable accession format, return the second part
        return parts[1]
    
    # Handle simple format with version number (e.g., CP101522.1)
    if '.' in subject_id and not ' ' in subject_id:
        return subject_id.split(':')[0]
    
    # Return the original ID if no known format is detected
    return subject_id

def extract_genome_info(subject_id):
    """
    Maps the subject ID to a simplified genome name using the global mapping.
    
    Args:
        subject_id (str): Subject ID from BLAST results
        
    Returns:
        str: Simplified genome name (e.g., "Genome_1")
    """
    global genome_mapping
    
    # Check if the mapping exists and contains this subject ID
    if genome_mapping and subject_id in genome_mapping:
        return genome_mapping[subject_id]
    
    # Fallback to original extraction method if no mapping exists
    return extract_genome_info_original(subject_id)

def create_heatmaps(blast_df, output_dir, prefix="adapted_", distance_metric="euclidean", linkage_method="average"):
    """
    Create two heatmap visualizations from the parsed BLAST data with hierarchical clustering:
    1. A presence/absence heatmap with dendrograms
    2. A percent identity heatmap with dendrograms
    
    Args:
        blast_df (pandas.DataFrame): DataFrame containing parsed BLAST data
        output_dir (str): Directory to save the output visualizations
        prefix (str): Prefix to add to output filenames
        distance_metric (str): Distance metric for clustering (euclidean, correlation, jaccard, etc.)
        linkage_method (str): Linkage method for hierarchical clustering (single, complete, average, ward, etc.)
    """
    global genome_mapping
    
    # Create genome name mapping
    genome_mapping = create_genome_name_mapping(blast_df)
    
    # Save the genome mapping to a separate file
    mapping_path = os.path.join(output_dir, f'{prefix}genome_name_mapping.csv')
    with open(mapping_path, 'w') as f:
        f.write("Original_ID,Cleaned_ID,Simple_Name\n")
        for original_id, simple_name in genome_mapping.items():
            cleaned_id = extract_genome_info_original(original_id)
            f.write(f"{original_id},{cleaned_id},{simple_name}\n")
            
    print(f"Genome name mapping saved to {mapping_path}")
    # Extract region and genome information
    blast_df['region'] = blast_df['query_id'].apply(extract_region_info)
    blast_df['genome'] = blast_df['subject_id'].apply(extract_genome_info)
    
    # Create a pivot table with regions as columns, genomes as rows, and percent identity as values
    pivot_data = blast_df.pivot_table(
        index='genome',  # organisms on y-axis
        columns='region',  # regions on x-axis
        values='percent_identity',
        aggfunc='mean'  # Using mean for regions/genomes with multiple hits
    )
    
    # Fill NaN values with 0 (no alignment)
    pivot_data = pivot_data.fillna(0)
    
    # Create presence/absence binary matrix
    presence_absence = pivot_data.copy()
    presence_absence[presence_absence > 0] = 1
    
    # Determine figure dimensions based on data size
    n_regions = len(pivot_data.columns)
    n_genomes = len(pivot_data.index)
    fig_width = max(16, n_regions * 0.4)  # Dynamic width based on number of regions
    fig_height = max(10, n_genomes * 0.4)  # Dynamic height based on number of genomes
    
    # 1. Create the presence/absence heatmap with hierarchical clustering and dendrograms
    
    # Perform hierarchical clustering on rows (genomes)
    row_linkage = hierarchy.linkage(
        pdist(presence_absence, metric=distance_metric),
        method=linkage_method
    )
    
    # Perform hierarchical clustering on columns (regions)
    col_linkage = hierarchy.linkage(
        pdist(presence_absence.T, metric=distance_metric),
        method=linkage_method
    )
    
    # Create a figure with gridspec for main heatmap and dendrograms
    fig = plt.figure(figsize=(fig_width + 3, fig_height + 3))
    gs = gridspec.GridSpec(2, 2, width_ratios=[0.15, 0.85], height_ratios=[0.15, 0.85])
    
    # Dendrogram for columns (regions) - top
    ax_col_dendrogram = fig.add_subplot(gs[0, 1])
    col_dendrogram = hierarchy.dendrogram(
        col_linkage, 
        color_threshold=0,
        ax=ax_col_dendrogram
    )
    ax_col_dendrogram.set_xticks([])
    ax_col_dendrogram.set_yticks([])
    ax_col_dendrogram.set_title('Region Clustering', fontsize=12)
    
    # Dendrogram for rows (genomes) - left
    ax_row_dendrogram = fig.add_subplot(gs[1, 0])
    row_dendrogram = hierarchy.dendrogram(
        row_linkage, 
        color_threshold=0,
        orientation='left',
        ax=ax_row_dendrogram
    )
    ax_row_dendrogram.set_xticks([])
    ax_row_dendrogram.set_yticks([])
    ax_row_dendrogram.set_title('Genome Clustering', fontsize=12, rotation=90, y=0.5, va='center')
    
    # Reorder data based on clustering
    row_order = row_dendrogram['leaves']
    col_order = col_dendrogram['leaves']
    presence_absence_clustered = presence_absence.iloc[row_order, col_order]
    
    # Main heatmap with clustered data
    ax_heatmap = fig.add_subplot(gs[1, 1])
    heatmap_pa = sns.heatmap(
        presence_absence_clustered,
        cmap=['white', 'darkblue'],  # Binary color scheme
        annot=False,     
        linewidths=0.5,
        cbar_kws={'label': 'Presence (1) / Absence (0)'},
        ax=ax_heatmap
    )
    
    # Improve axis labels
    plt.suptitle('Presence/Absence of Divergent Regions Across Organisms', fontsize=16, y=0.98)
    ax_heatmap.set_xlabel('Divergent Regions (Coordinates)', fontsize=14)
    ax_heatmap.set_ylabel('Organisms', fontsize=14)
    
    # Rotate x-axis labels for better readability
    plt.setp(ax_heatmap.get_xticklabels(), rotation=90, ha='center', fontsize=9)
    plt.setp(ax_heatmap.get_yticklabels(), fontsize=10)
    
    # Adjust layout
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    
    # Save the figure
    presence_absence_path = os.path.join(output_dir, f'{prefix}blast_presence_absence_heatmap.png')
    plt.savefig(presence_absence_path, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"Presence/Absence heatmap saved to {presence_absence_path}")
    
    # 2. Create the percent identity heatmap with hierarchical clustering and dendrograms
    
    # Create a mask for cells with 0% identity
    mask = pivot_data == 0
    
    # Perform hierarchical clustering on rows (genomes)
    # For percent identity, we might want to use correlation or cosine distance
    row_linkage = hierarchy.linkage(
        pdist(pivot_data.fillna(0), metric=distance_metric),
        method=linkage_method
    )
    
    # Perform hierarchical clustering on columns (regions)
    col_linkage = hierarchy.linkage(
        pdist(pivot_data.fillna(0).T, metric=distance_metric),
        method=linkage_method
    )
    
    # Create a figure with gridspec for main heatmap and dendrograms
    fig = plt.figure(figsize=(fig_width + 3, fig_height + 3))
    gs = gridspec.GridSpec(2, 2, width_ratios=[0.15, 0.85], height_ratios=[0.15, 0.85])
    
    # Dendrogram for columns (regions) - top
    ax_col_dendrogram = fig.add_subplot(gs[0, 1])
    col_dendrogram = hierarchy.dendrogram(
        col_linkage, 
        color_threshold=0,
        ax=ax_col_dendrogram
    )
    ax_col_dendrogram.set_xticks([])
    ax_col_dendrogram.set_yticks([])
    ax_col_dendrogram.set_title('Region Clustering', fontsize=12)
    
    # Dendrogram for rows (genomes) - left
    ax_row_dendrogram = fig.add_subplot(gs[1, 0])
    row_dendrogram = hierarchy.dendrogram(
        row_linkage, 
        color_threshold=0,
        orientation='left',
        ax=ax_row_dendrogram
    )
    ax_row_dendrogram.set_xticks([])
    ax_row_dendrogram.set_yticks([])
    ax_row_dendrogram.set_title('Genome Clustering', fontsize=12, rotation=90, y=0.5, va='center')
    
    # Reorder data based on clustering
    row_order = row_dendrogram['leaves']
    col_order = col_dendrogram['leaves']
    pivot_data_clustered = pivot_data.iloc[row_order, col_order]
    mask_clustered = mask.iloc[row_order, col_order]
    
    # Main heatmap with clustered data
    ax_heatmap = fig.add_subplot(gs[1, 1])
    heatmap_id = sns.heatmap(
        pivot_data_clustered,
        cmap='Blues',  # Blue color scale
        annot=True,    # Show percent identity values in cells
        fmt=".1f",     # Format percent identity to 1 decimal place
        annot_kws={"size": 7, "color": "black"},  # Adjust annotation appearance
        linewidths=0.5,
        cbar_kws={'label': 'Percent Identity (%)'},
        mask=mask_clustered,  # Hide cells with 0% identity
        ax=ax_heatmap
    )
    
    # Function to show only non-zero values
    for text in heatmap_id.texts:
        if text.get_text():  # Check if text is not empty
            value = float(text.get_text())
            if value == 0:
                text.set_text("")  # Remove text for cells with 0% identity
    
    # Improve axis labels
    plt.suptitle('Percent Identity of Divergent Regions Across Organisms', fontsize=16, y=0.98)
    ax_heatmap.set_xlabel('Divergent Regions (Coordinates)', fontsize=14)
    ax_heatmap.set_ylabel('Organisms', fontsize=14)
    
    # Rotate x-axis labels for better readability
    plt.setp(ax_heatmap.get_xticklabels(), rotation=90, ha='center', fontsize=9)
    plt.setp(ax_heatmap.get_yticklabels(), fontsize=10)
    
    # Adjust layout
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    
    # Save the figure
    percent_identity_path = os.path.join(output_dir, f'{prefix}blast_percent_identity_heatmap.png')
    plt.savefig(percent_identity_path, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"Percent identity heatmap saved to {percent_identity_path}")
    
    # Save the processed data for future reference
    data_path = os.path.join(output_dir, f'{prefix}blast_heatmap_data.csv')
    pivot_data.to_csv(data_path)
    print(f"Processed data saved to {data_path}")

# Global variable to store genome mapping
genome_mapping = {}

def main():
    """Main function to execute the heatmap generation workflow."""
    # Define file paths
    blast_results_path = 'results/blast_results.txt'
    output_dir = 'results'
    
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    # Parse the BLAST results
    print(f"Parsing BLAST results from {blast_results_path}...")
    blast_df = parse_blast_results(blast_results_path)
    
    if blast_df is None or blast_df.empty:
        print("Error: Unable to parse BLAST results or file is empty.")
        return
    
    print(f"Found {len(blast_df)} BLAST alignment records.")
    
    # Create and save the heatmaps
    print("Generating heatmap visualizations with optimized identifier extraction...")
    create_heatmaps(blast_df, output_dir)
    
    print("Process completed successfully.")
    print("Note: This script is specifically adapted to handle identifiers in the format NZ_CP101522.1:XXXXX-XXXXX.")
    print("Genome identifiers have been simplified to Genome_1, Genome_2, etc.")
    print("Zero percent identity cells have been masked out from the heatmap.")
    print("Genome name mapping has been saved for reference.")

if __name__ == "__main__":
    main()

