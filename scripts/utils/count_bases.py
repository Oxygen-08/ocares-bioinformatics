def count_bases(fasta_file):
    """Counts the total number of bases in a FASTA file."""
    total_bases = 0
    with open(fasta_file, 'r') as file:
        for line in file:
            if not line.startswith('>'):  # Skip header lines
                total_bases += len(line.strip())  # Add length of sequence line
    return total_bases

# Example usage
fasta_file = 'genome.fasta'  # Replace with your file path
print(f"Total number of bases: {count_bases(fasta_file)}")
