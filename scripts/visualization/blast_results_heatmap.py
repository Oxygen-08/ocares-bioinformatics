#!/usr/bin/env python3
"""
Generate heatmap visualizations from BLAST results data.
This script reads BLAST results, extracts query regions, subject genomes, 
and percent identity to create two heatmaps:
1. A presence/absence heatmap (binary: 0/1)
2. A percent identity heatmap using blue color scale
"""

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path

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
    Extract genomic region information from the query ID.
    Adjust this function based on your specific query ID format.
    
    Args:
        query_id (str): Query ID from BLAST results
        
    Returns:
        str: Extracted region information
    """
    # Example implementation - modify according to your query ID format
    # This assumes query IDs contain region information (like chromosome or contig)
    parts = str(query_id).split('|')
    if len(parts) > 1:
        return parts[0]
    
    # Alternative parsing if query IDs are in a different format
    parts = str(query_id).split('_')
    if len(parts) > 0:
        return parts[0]
    
    return str(query_id)

def extract_genome_info(subject_id):
    """
    Extract genome identifier from the subject ID.
    Adjust this function based on your specific subject ID format.
    
    Args:
        subject_id (str): Subject ID from BLAST results
        
    Returns:
        str: Extracted genome information
    """
    # Example implementation - modify according to your subject ID format
    # This assumes subject IDs contain genome identifiers
    parts = str(subject_id).split('|')
    if len(parts) > 1:
        return parts[0]
        
    # Alternative parsing if subject IDs are in a different format
    parts = str(subject_id).split('_')
    if len(parts) > 0:
        return parts[0]
    
    return str(subject_id)

def create_heatmaps(blast_df, output_dir):
    """
    Create two heatmap visualizations from the parsed BLAST data:
    1. A presence/absence heatmap
    2. A percent identity heatmap
    
    Args:
        blast_df (pandas.DataFrame): DataFrame containing parsed BLAST data
        output_dir (str): Directory to save the output visualizations
    """
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
    
    # 1. Create the presence/absence heatmap
    plt.figure(figsize=(14, 10))
    
    # Use a categorical color map (white for absence, dark blue for presence)
    heatmap_pa = sns.heatmap(
        presence_absence,
        cmap=['white', 'darkblue'],  # Binary color scheme
        annot=False,     
        linewidths=0.5,
        cbar_kws={'label': 'Presence (1) / Absence (0)'}
    )
    
    # Improve axis labels
    plt.title('Presence/Absence of Divergent Regions Across Organisms', fontsize=16, pad=20)
    plt.xlabel('Divergent Regions', fontsize=14, labelpad=15)
    plt.ylabel('Organisms', fontsize=14, labelpad=15)
    
    # Rotate x-axis labels for better readability
    plt.xticks(rotation=45, ha='right', fontsize=10)
    plt.yticks(fontsize=10)
    
    # Adjust layout to prevent clipping of labels
    plt.tight_layout()
    
    # Save the figure
    presence_absence_path = os.path.join(output_dir, 'blast_presence_absence_heatmap.png')
    plt.savefig(presence_absence_path, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"Presence/Absence heatmap saved to {presence_absence_path}")
    
    # 2. Create the percent identity heatmap
    plt.figure(figsize=(14, 10))
    
    # Use a blue color scale for percent identity
    heatmap_id = sns.heatmap(
        pivot_data,
        cmap='Blues',  # Blue color scale
        annot=False,
        linewidths=0.5,
        cbar_kws={'label': 'Percent Identity (%)'}
    )
    
    # Improve axis labels
    plt.title('Percent Identity of Divergent Regions Across Organisms', fontsize=16, pad=20)
    plt.xlabel('Divergent Regions', fontsize=14, labelpad=15)
    plt.ylabel('Organisms', fontsize=14, labelpad=15)
    
    # Rotate x-axis labels for better readability
    plt.xticks(rotation=45, ha='right', fontsize=10)
    plt.yticks(fontsize=10)
    
    # Adjust layout to prevent clipping of labels
    plt.tight_layout()
    
    # Save the figure
    percent_identity_path = os.path.join(output_dir, 'blast_percent_identity_heatmap.png')
    plt.savefig(percent_identity_path, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"Percent identity heatmap saved to {percent_identity_path}")

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
    print("Generating heatmap visualizations...")
    create_heatmaps(blast_df, output_dir)
    
    print("Process completed successfully.")

if __name__ == "__main__":
    main()

