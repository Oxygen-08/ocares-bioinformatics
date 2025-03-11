#!/usr/bin/env python3
"""
BLAST Statistics Visualization Tool

This script generates statistical visualizations for BLAST results including:
1. Distribution plots for percent identity
2. Distribution of alignment lengths
3. E-value distribution (log-transformed)
4. Box plots for percent identity across genomes
5. Summary statistics table

All plots are saved in the 'results/statistics' directory.
"""

import os
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.gridspec import GridSpec
import math
import warnings
warnings.filterwarnings('ignore')

# Set seaborn style for attractive plots
sns.set(style="whitegrid")
plt.rcParams['figure.figsize'] = [10, 6]
plt.rcParams['figure.dpi'] = 300
sns.set_palette("viridis")

def ensure_dir(directory):
    """Create directory if it doesn't exist."""
    if not os.path.exists(directory):
        os.makedirs(directory)

def read_blast_results(blast_file):
    """
    Read BLAST results file and return a properly formatted DataFrame.
    
    Assumes standard BLAST output format 6 (tabular):
    qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore
    """
    columns = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 
               'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore']
    
    try:
        df = pd.read_csv(blast_file, sep='\t', header=None, names=columns)
        print(f"Successfully read {len(df)} BLAST alignments")
        return df
    except Exception as e:
        print(f"Error reading BLAST file: {e}")
        sys.exit(1)

def extract_genome_info(sseqid):
    """Extract a clean genome identifier from the subject sequence ID."""
    # Handle different formats of subject IDs
    if '|' in sseqid:
        # Format like gi|123456789|ref|NC_123456.1|
        parts = sseqid.split('|')
        for part in parts:
            if part.startswith(('NC_', 'NZ_', 'NG_')):
                return part
        return parts[3] if len(parts) > 3 else parts[1]
    else:
        # Format like NC_123456.1 or NZ_CP101522.1
        return sseqid.split()[0].split(':')[0]

def plot_percent_identity_distribution(df, output_dir):
    """Generate percent identity distribution plot."""
    plt.figure(figsize=(10, 6))
    
    ax = sns.histplot(df['pident'], bins=30, kde=True, color='mediumseagreen')
    ax.set_title('Distribution of Percent Identity', fontsize=16)
    ax.set_xlabel('Percent Identity', fontsize=14)
    ax.set_ylabel('Frequency', fontsize=14)
    
    # Add vertical lines for mean and median
    mean_pident = df['pident'].mean()
    median_pident = df['pident'].median()
    
    plt.axvline(mean_pident, color='darkblue', linestyle='--', 
                label=f'Mean: {mean_pident:.2f}%')
    plt.axvline(median_pident, color='crimson', linestyle='-', 
                label=f'Median: {median_pident:.2f}%')
    
    plt.legend(fontsize=12)
    
    # Add text annotations for percentages
    high_identity = (df['pident'] >= 95).mean() * 100
    perfect_identity = (df['pident'] == 100).mean() * 100
    
    plt.annotate(f"{high_identity:.1f}% of alignments have ≥95% identity", 
                xy=(95, plt.ylim()[1] * 0.9), 
                xytext=(85, plt.ylim()[1] * 0.95),
                arrowprops=dict(arrowstyle="->", color='black'))
    
    plt.annotate(f"{perfect_identity:.1f}% of alignments have 100% identity", 
                xy=(100, plt.ylim()[1] * 0.7), 
                xytext=(90, plt.ylim()[1] * 0.75),
                arrowprops=dict(arrowstyle="->", color='black'))
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'percent_identity_distribution.png'))
    plt.close()

def plot_alignment_length_distribution(df, output_dir):
    """Generate alignment length distribution plot."""
    plt.figure(figsize=(10, 6))
    
    # Create log-scale for x-axis if lengths vary significantly
    if df['length'].max() / df['length'].min() > 100:
        plt.xscale('log')
        bins = np.logspace(np.log10(df['length'].min()), np.log10(df['length'].max()), 30)
    else:
        bins = 30
    
    ax = sns.histplot(df['length'], bins=bins, kde=True, color='darkorange')
    ax.set_title('Distribution of Alignment Lengths', fontsize=16)
    ax.set_xlabel('Alignment Length (bp)', fontsize=14)
    ax.set_ylabel('Frequency', fontsize=14)
    
    # Add vertical lines and text for statistics
    mean_len = df['length'].mean()
    median_len = df['length'].median()
    
    plt.axvline(mean_len, color='darkblue', linestyle='--', 
                label=f'Mean: {mean_len:.1f} bp')
    plt.axvline(median_len, color='crimson', linestyle='-', 
                label=f'Median: {median_len:.1f} bp')
    
    # Annotate proportion of short alignments
    short_alignments = (df['length'] < 50).mean() * 100
    
    if short_alignments > 0:
        plt.annotate(f"{short_alignments:.1f}% of alignments are <50 bp", 
                    xy=(50, plt.ylim()[1] * 0.8), 
                    xytext=(max(100, df['length'].median()/2), plt.ylim()[1] * 0.85),
                    arrowprops=dict(arrowstyle="->", color='black'))
    
    plt.legend(fontsize=12)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'alignment_length_distribution.png'))
    plt.close()

def plot_evalue_distribution(df, output_dir):
    """Generate E-value distribution plot (log-transformed)."""
    plt.figure(figsize=(10, 6))
    
    # Handle E-values of 0 by setting them to a very small number
    df_plot = df.copy()
    df_plot.loc[df_plot['evalue'] == 0, 'evalue'] = 1e-300
    
    # Log-transform E-values
    log_evalues = -np.log10(df_plot['evalue'])
    
    ax = sns.histplot(log_evalues, bins=30, kde=True, color='royalblue')
    ax.set_title('Distribution of E-values (Log-transformed)', fontsize=16)
    ax.set_xlabel('-log10(E-value)', fontsize=14)
    ax.set_ylabel('Frequency', fontsize=14)
    
    # Add vertical lines for common significance thresholds
    plt.axvline(-np.log10(1e-50), color='red', linestyle='--', 
                label='E-value = 1e-50')
    plt.axvline(-np.log10(1e-10), color='orange', linestyle='--', 
                label='E-value = 1e-10')
    
    # Calculate percentage of zero E-values
    zero_evalues = (df['evalue'] == 0).mean() * 100
    
    # Annotate the plot
    plt.annotate(f"{zero_evalues:.1f}% of alignments have E-value = 0", 
                xy=(plt.xlim()[1] * 0.9, plt.ylim()[1] * 0.9), 
                xytext=(plt.xlim()[1] * 0.6, plt.ylim()[1] * 0.95),
                arrowprops=dict(arrowstyle="->", color='black'))
    
    plt.legend(fontsize=12)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'evalue_distribution.png'))
    plt.close()

def plot_identity_by_genome(df, output_dir):
    """Generate box plots for percent identity across genomes."""
    # Extract genome information from subject IDs
    df['genome'] = df['sseqid'].apply(extract_genome_info)
    
    # Count alignments per genome for sorting
    genome_counts = df['genome'].value_counts()
    
    # Get top genomes (to avoid overcrowding the plot)
    top_genomes = genome_counts.head(15).index.tolist()
    
    # Filter data for top genomes
    df_plot = df[df['genome'].isin(top_genomes)].copy()
    
    # Sort by median percent identity
    genome_median = df_plot.groupby('genome')['pident'].median().sort_values()
    ordered_genomes = genome_median.index.tolist()
    
    plt.figure(figsize=(12, 8))
    ax = sns.boxplot(x='genome', y='pident', data=df_plot, 
                    order=ordered_genomes, palette='viridis')
    
    ax.set_title('Percent Identity by Genome', fontsize=16)
    ax.set_xlabel('Genome', fontsize=14)
    ax.set_ylabel('Percent Identity', fontsize=14)
    
    # Rotate x-axis labels for better readability
    plt.xticks(rotation=45, ha='right')
    
    # Add count information
    for i, genome in enumerate(ordered_genomes):
        count = genome_counts[genome]
        ax.annotate(f"n={count}", 
                   xy=(i, df_plot[df_plot['genome'] == genome]['pident'].min() - 1),
                   ha='center')
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'identity_by_genome.png'))
    plt.close()

def create_summary_table(df, output_dir):
    """Generate summary statistics table as a figure."""
    # Calculate overall statistics
    stats = {
        "Total Alignments": len(df),
        "Unique Query Sequences": df['qseqid'].nunique(),
        "Unique Subject Sequences": df['sseqid'].nunique(),
        "Mean Percent Identity": f"{df['pident'].mean():.2f}%",
        "Median Percent Identity": f"{df['pident'].median():.2f}%",
        "Mean Alignment Length": f"{df['length'].mean():.1f} bp",
        "Median Alignment Length": f"{df['length'].median():.1f} bp",
        "Alignments with E-value = 0": f"{(df['evalue'] == 0).sum()} ({(df['evalue'] == 0).mean()*100:.1f}%)",
        "Alignments with ≥95% Identity": f"{(df['pident'] >= 95).sum()} ({(df['pident'] >= 95).mean()*100:.1f}%)",
        "Alignments with 100% Identity": f"{(df['pident'] == 100).sum()} ({(df['pident'] == 100).mean()*100:.1f}%)",
        "Short Alignments (<50 bp)": f"{(df['length'] < 50).sum()} ({(df['length'] < 50).mean()*100:.1f}%)"
    }
    
    # Create a DataFrame for easier plotting
    stats_df = pd.DataFrame(list(stats.items()), columns=['Metric', 'Value'])
    
    # Create figure and plot
    fig, ax = plt.figure(figsize=(10, 6)), plt.subplot(111)
    ax.axis('off')
    
    # Hide axes
    ax.axis('tight')
    ax.axis('off')
    
    # Create table
    table = ax.table(cellText=stats_df.values, colLabels=stats_df.columns,
                    loc='center', cellLoc='left', colWidths=[0.6, 0.4])
    
    # Style the table
    table.auto_set_font_size(False)
    table.set_fontsize(12)
    table.scale(1.2, 1.8)
    
    # Set title
    plt.title("BLAST Results Summary Statistics", fontsize=16, y=0.8)
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'summary_statistics.png'))
    plt.close()

def create_combined_dashboard(df, output_dir):
    """Create a combined dashboard with all visualizations."""
    # Create a large figure
    fig = plt.figure(figsize=(20, 16))
    gs = GridSpec(2, 2, figure=fig)
    
    # Plot percent identity distribution
    ax1 = fig.add_subplot(gs[0, 0])
    df_plot = df.copy()
    sns.histplot(df_plot['pident'], bins=30, kde=True, color='mediumseagreen', ax=ax1)
    ax1.set_title('Distribution of Percent Identity', fontsize=14)
    ax1.set_xlabel('Percent Identity', fontsize=12)
    ax1.set_ylabel('Frequency', fontsize=12)
    
    # Plot alignment length distribution
    ax2 = fig.add_subplot(gs[0, 1])
    if df['length'].max() / df['length'].min() > 100:
        ax2.set_xscale('log')
        bins = np.logspace(np.log10(df['length'].min()), np.log10(df['length'].max()), 30)
    else:
        bins = 30
    sns.histplot(df['length'], bins=bins, kde=True, color='darkorange', ax=ax2)
    ax2.set_title('Distribution of Alignment Lengths', fontsize=14)
    ax2.set_xlabel('Alignment Length (bp)', fontsize=12)
    ax2.set_ylabel('Frequency', fontsize=12)
    
    # Plot E-value distribution
    ax3 = fig.add_subplot(gs[1, 0])
    df_plot = df.copy()
    df_plot.loc[df_plot['evalue'] == 0, 'evalue'] = 1e-300
    log_evalues

