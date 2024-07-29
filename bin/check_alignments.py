#!/usr/bin/env python
import argparse
import numpy as np

def read_fasta(file_path):
    """Read a FASTA file and return a dictionary of sequences."""
    sequences = {}
    with open(file_path, 'r') as file:
        name, seq = None, []
        for line in file:
            line = line.strip()
            if line.startswith(">"):
                if name:
                    sequences[name] = ''.join(seq)
                name, seq = line[1:], []
            else:
                seq.append(line)
        if name:
            sequences[name] = ''.join(seq)
    return sequences

def write_fasta(sequences, file_path):
    """Write sequences to a FASTA file."""
    with open(file_path, 'w') as file:
        for name, seq in sequences.items():
            file.write(f">{name}\n")
            file.write(f"{seq}\n")

def remove_monomorphic_columns(seq_array):
    """Remove columns that are identical across all sequences."""
    return seq_array[:, ~(seq_array == seq_array[0]).all(axis=0)]

def nucleotides_to_integers(seq_array):
    """Convert nucleotide sequences to integer arrays."""
    char_to_int = {char: idx for idx, char in enumerate('ACGT')}
    return np.vectorize(char_to_int.get)(seq_array)

def hamming_distance_vectorized(seq_array):
    """Calculate the Hamming distances between all pairs of sequences using optimized vectorized operations."""
    n = seq_array.shape[0]
    diff_matrix = seq_array[:, np.newaxis, :] != seq_array[np.newaxis, :, :]
    dist_matrix = np.sum(diff_matrix, axis=2)
    
    # Get the upper triangle of the distance matrix as a 1D array
    distances = dist_matrix[np.triu_indices(n, k=1)]
    return distances

def analyze_sequences(fasta_file, num_tips, prefix):
    # Read sequences from FASTA file
    sequences = read_fasta(fasta_file)
    sequence_list = list(sequences.values())

    # Convert sequences to numpy array for vectorized operations
    seq_array = np.array([list(seq) for seq in sequence_list])

    # Remove monomorphic columns
    seq_array = remove_monomorphic_columns(seq_array)

    # Convert nucleotides to integers
    seq_array = nucleotides_to_integers(seq_array)

    # Calculate number of unique sequences and duplicates
    unique_seqs, counts = np.unique(seq_array, axis=0, return_counts=True)
    num_unique = len(unique_seqs)
    num_duplicates = len(sequence_list) - num_unique

    # Calculate Hamming distances between sequences
    distances = hamming_distance_vectorized(unique_seqs)
    min_distance = np.min(distances)
    max_distance = np.max(distances)
    avg_distance = np.mean(distances)

    # Write the metrics to a file
    report_file = f"{prefix}_report.txt"
    with open(report_file, 'w') as f:
        f.write(f"Number of unique sequences: {num_unique}\n")
        f.write(f"Number of duplicate sequences: {num_duplicates}\n")
        f.write(f"Minimum Hamming distance: {min_distance}\n")
        f.write(f"Maximum Hamming distance: {max_distance}\n")
        f.write(f"Average Hamming distance: {avg_distance}\n")

def main():
    parser = argparse.ArgumentParser(description="Analyze alignments and gather sequence metrics.")
    parser.add_argument('--fasta', type=str, required=True, help="Input FASTA file")
    parser.add_argument('--tips', type=int, required=True, help="Expected number of tips")
    parser.add_argument('--prefix', type=str, required=True, help="Output file prefix")
    args = parser.parse_args()

    analyze_sequences(args.fasta, args.tips, args.prefix)

if __name__ == "__main__":
    main()
