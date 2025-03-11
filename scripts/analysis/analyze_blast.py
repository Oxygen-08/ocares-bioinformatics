#!/usr/bin/env python3
import os
import sys
import pandas as pd
import numpy as np
import re

def read_blast_file(file_path):
    """
    Read BLAST output file into a pandas DataFrame.
    
    Args:
        file_path (str): Path to the BLAST output file
        
    Returns:
        pandas.DataFrame: DataFrame containing BLAST results or None if an error occurs
    """
    try:
        if not os.path.exists(file_path):
            print(f"Error: File '{file_path}' not found.")
            return None
            
        print(f"Reading BLAST output from '{file_path}'...")
        
        # Define column names based on standard BLAST tabular output
        columns = [
            'query_id', 'subject_id', 'percent_identity', 'alignment_length',
            'mismatches', 'gap_opens', 'q_start', 'q_end', 's_start', 's_end',
            'e_value', 'bit_score'
        ]
        
        # Try to read the file with different delimiters
        try:
            # First try tab-delimited (standard BLAST output)
            df = pd.read_csv(file_path, sep='\t', header=None, names=columns)
        except pd.errors.ParserError:
            # If that fails, try whitespace-delimited
            df = pd.read_csv(file_path, sep='\s+', header=None, names=columns)
            
        # Check if the DataFrame is empty
        if df.empty:
            print("Warning: The BLAST output file is empty.")
            return None
            
        print(f"Successfully loaded {len(df)} BLAST result entries.")
        return df
        
    except Exception as e:
        print(f"Error reading BLAST file: {str(e)}")
        return None

def clean_data(df):
    """
    Clean and preprocess the BLAST data.
    
    Args:
        df (pandas.DataFrame): DataFrame containing BLAST results
        
    Returns:
        pandas.DataFrame: Cleaned DataFrame
    """
    if df is None:
        return None
        
    try:
        # Make a copy to avoid modifying the original
        clean_df = df.copy()
        
        # Handle non-numeric percent identity values (like asterisks)
        if 'percent_identity' in clean_df.columns:
            # Replace non-numeric values with NaN
            clean_df['percent_identity'] = pd.to_numeric(
                clean_df['percent_identity'], errors='coerce')
            
            # Count how many were converted to NaN
            nan_count = clean_df['percent_identity'].isna().sum()
            if nan_count > 0:
                print(f"Warning: {nan_count} entries had non-numeric percent identity values.")
        
        # Remove rows with missing critical values
        original_rows = len(clean_df)
        clean_df = clean_df.dropna(subset=['query_id', 'subject_id'])
        
        if len(clean_df) < original_rows:
            print(f"Removed {original_rows - len(clean_df)} rows with missing query or subject IDs.")
        
        # Convert numeric columns to appropriate types
        numeric_cols = ['alignment_length', 'mismatches', 'gap_opens', 
                        'q_start', 'q_end', 's_start', 's_end', 'bit_score']
                        
        for col in numeric_cols:
            if col in clean_df.columns:
                clean_df[col] = pd.to_numeric(clean_df[col], errors='coerce')
        
        # Handle e-value scientific notation
        if 'e_value' in clean_df.columns:
            clean_df['e_value'] = pd.to_numeric(clean_df['e_value'], errors='coerce')
            
        return clean_df
        
    except Exception as e:
        print(f"Error cleaning data: {str(e)}")
        return df  # Return original if cleaning fails

def analyze_blast_data(df):
    """
    Analyze BLAST data and print statistics and potential issues.
    
    Args:
        df (pandas.DataFrame): DataFrame containing BLAST results
    """
    if df is None:
        return
        
    print("\n" + "="*80)
    print("BLAST DATA ANALYSIS SUMMARY")
    print("="*80)
    
    # 1. Basic dataset information
    print(f"\nTotal number of alignments: {len(df):,}")
    
    # 2. Query and subject statistics
    query_count = df['query_id'].nunique()
    subject_count = df['subject_id'].nunique()
    print(f"Number of unique query sequences: {query_count:,}")
    print(f"Number of unique subject sequences: {subject_count:,}")
    
    # 3. Percent identity statistics
    if 'percent_identity' in df.columns:
        id_stats = df['percent_identity'].describe()
        print("\nPercent Identity Statistics:")
        print(f"  Min: {id_stats['min']:.2f}%")
        print(f"  Max: {id_stats['max']:.2f}%")
        print(f"  Mean: {id_stats['mean']:.2f}%")
        print(f"  Median: {id_stats['50%']:.2f}%")
        
        # Count alignments in different identity ranges
        id_ranges = [
            (0, 25), (25, 50), (50, 75), (75, 90), 
            (90, 95), (95, 99), (99, 100), (100, 100.1)
        ]
        
        print("\nPercent Identity Distribution:")
        for low, high in id_ranges:
            if high == 100.1:  # Handle exact 100% case
                count = (df['percent_identity'] == 100).sum()
                print(f"  Exactly 100%: {count:,} alignments ({count/len(df)*100:.2f}%)")
            else:
                count = ((df['percent_identity'] >= low) & (df['percent_identity'] < high)).sum()
                print(f"  {low}% to <{high}%: {count:,} alignments ({count/len(df)*100:.2f}%)")
    
    # 4. E-value statistics
    if 'e_value' in df.columns:
        print("\nE-value Statistics:")
        e_value_zero = (df['e_value'] == 0).sum()
        print(f"  Alignments with E-value = 0: {e_value_zero:,} ({e_value_zero/len(df)*100:.2f}%)")
        
        # Count alignments in different E-value ranges (using log scale)
        if df['e_value'].min() > 0:
            print(f"  Minimum non-zero E-value: {df['e_value'].min():.2e}")
        
        significant_hits = (df['e_value'] <= 1e-50).sum()
        print(f"  Highly significant hits (E ≤ 1e-50): {significant_hits:,} ({significant_hits/len(df)*100:.2f}%)")
    
    # 5. Alignment length statistics
    if 'alignment_length' in df.columns:
        length_stats = df['alignment_length'].describe()
        print("\nAlignment Length Statistics:")
        print(f"  Min: {length_stats['min']:.0f} bp")
        print(f"  Max: {length_stats['max']:.0f} bp")
        print(f"  Mean: {length_stats['mean']:.2f} bp")
        print(f"  Median: {length_stats['50%']:.0f} bp")
    
    # 6. Identify potential issues
    print("\n" + "="*80)
    print("POTENTIAL DATA ISSUES")
    print("="*80)
    
    # Check for missing values in each column
    missing_values = df.isna().sum()
    if missing_values.sum() > 0:
        print("\nMissing Values:")
        for col, count in missing_values.items():
            if count > 0:
                print(f"  {col}: {count:,} missing values ({count/len(df)*100:.2f}%)")
    
    # Check for unusually low percent identity values
    if 'percent_identity' in df.columns:
        low_identity = (df['percent_identity'] < 30).sum()
        if low_identity > 0:
            print(f"\nLow Percent Identity:")
            print(f"  {low_identity:,} alignments have percent identity < 30% ({low_identity/len(df)*100:.2f}%)")
            print("  This may indicate spurious matches or divergent sequences.")
    
    # Check for very short alignments
    if 'alignment_length' in df.columns:
        short_alignments = (df['alignment_length'] < 50).sum()
        if short_alignments > 0:
            print(f"\nShort Alignments:")
            print(f"  {short_alignments:,} alignments are shorter than 50 bp ({short_alignments/len(df)*100:.2f}%)")
            print("  Very short alignments may be spurious or not biologically meaningful.")
    
    # Check for unusual gap patterns
    if 'gap_opens' in df.columns and 'alignment_length' in df.columns:
        high_gap_ratio = ((df['gap_opens'] / df['alignment_length']) > 0.1).sum()
        if high_gap_ratio > 0:
            print(f"\nHigh Gap Ratio:")
            print(f"  {high_gap_ratio:,} alignments have gaps in >10% of their length ({high_gap_ratio/len(df)*100:.2f}%)")
            print("  High gap frequency may indicate problematic alignments.")
    
    # Check for query ID and subject ID patterns
    print("\n" + "="*80)
    print("SEQUENCE IDENTIFIER ANALYSIS")
    print("="*80)
    
    # Analyze query ID patterns
    query_patterns = analyze_identifier_patterns(df['query_id'])
    print("\nQuery ID Patterns:")
    for pattern, count in query_patterns.items():
        print(f"  {pattern}: {count:,} sequences ({count/query_count*100:.2f}%)")
    
    # Analyze subject ID patterns
    subject_patterns = analyze_identifier_patterns(df['subject_id'])
    print("\nSubject ID Patterns:")
    for pattern, count in subject_patterns.items():
        print(f"  {pattern}: {count:,} sequences ({count/subject_count*100:.2f}%)")
    
    # Print some example IDs
    print("\nExample Query IDs:")
    for id in df['query_id'].drop_duplicates().head(5).tolist():
        print(f"  {id}")
    
    print("\nExample Subject IDs:")
    for id in df['subject_id'].drop_duplicates().head(5).tolist():
        print(f"  {id}")

def analyze_identifier_patterns(id_series):
    """
    Analyze patterns in sequence identifiers.
    
    Args:
        id_series (pandas.Series): Series containing sequence identifiers
        
    Returns:
        dict: Dictionary with pattern types and their counts
    """
    patterns = {}
    
    # Common identifier prefixes
    prefix_patterns = {
        'gi|': r'^gi\|',
        'NZ_': r'^NZ_',
        'NC_': r'^NC_',
        'AP_': r'^AP_',
        'CP_': r'^CP_',
        'NP_': r'^NP_',
        'WP_': r'^WP_',
        'YP_': r'^YP_',
        'XP_': r'^XP_',
        'sp|': r'^sp\|',
        'tr|': r'^tr\|',
        'lcl|': r'^lcl\|',
        'pdb|': r'^pdb\|',
        'ref|': r'^ref\|',
        'gb|': r'^gb\|',
        'emb|': r'^emb\|',
        'dbj|': r'^dbj\|',
    }
    
    unique_ids = id_series.drop_duplicates()
    
    # Count occurrences of each pattern
    for name, pattern in prefix_patterns.items():
        count = unique_ids.str.contains(pattern, regex=True).sum()
        if count > 0:
            patterns[name] = count
    
    # Count ids that don't match any known pattern
    all_patterns = '|'.join(prefix_patterns.values())
    unknown_count = (~unique_ids.str.contains(all_patterns, regex=True)).sum()
    if unknown_count > 0:
        patterns['unknown/other'] = unknown_count
    
    return patterns

def main():
    """
    Main function to run the BLAST analysis.
    """
    # Check if file path is provided as command-line argument
    if len(sys.argv) > 1:
        blast_file = sys.argv[1]
    else:
        # Default to the blast_results.txt file in the results directory
        blast_file = os.path.join('results', 'blast_results.txt')
    
    # Read and analyze the BLAST data
    df = read_blast_file(blast_file)
    if df is not None:
        clean_df = clean_data(df)
        if clean_df is not None:
            analyze_blast_data(clean_df)
            
            # Save the DataFrame to a CSV file for further analysis
            output_csv = os.path.splitext(blast_file)[0] + '_analyzed.csv'
            clean_df.to_csv(output_csv, index=False)
            print(f"\nCleaned data saved to: {output_csv}")
        else:
            print("Could not clean the BLAST data. Analysis aborted.")
    else:
        print("Could not read the BLAST data. Analysis aborted.")

if __name__ == "__main__":
    main()

