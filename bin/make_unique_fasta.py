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

def introduce_random_mutation(sequence):
    """Introduce a random mutation in the sequence."""
    pos = np.random.randint(len(sequence))
    bases = ['A', 'C', 'G', 'T']
    bases.remove(sequence[pos])
    new_base = np.random.choice(bases)
    return sequence[:pos] + new_base + sequence[pos + 1:]

def make_sequences_unique(seq_array):
    """Introduce random mutations to make sequences unique."""
    unique_seqs, counts = np.unique(seq_array, axis=0, return_counts=True)
    num_unique = len(unique_seqs)

    while num_unique < seq_array.shape[0]:
        for i, count in enumerate(counts):
            if count > 1:
                seq_array[i] = list(introduce_random_mutation(''.join(seq_array[i])))
        unique_seqs, counts = np.unique(seq_array, axis=0, return_counts=True)
        num_unique = len(unique_seqs)

    return seq_array

def ensure_unique_sequences(fasta_file, prefix):
    # Read sequences from FASTA file
    sequences = read_fasta(fasta_file)
    sequence_list = list(sequences.values())

    # Convert sequences to numpy array for vectorized operations
    seq_array = np.array([list(seq) for seq in sequence_list])

    # Ensure all sequences are unique
    seq_array = make_sequences_unique(seq_array)

    # Update the sequences with the forced unique ones
    for i, key in enumerate(sequences.keys()):
        sequences[key] = ''.join(seq_array[i])

    # Write the updated sequences to a new FASTA file
    validated_fasta_file = f"{prefix}.validated.fasta"
    write_fasta(sequences, validated_fasta_file)

def main():
    parser = argparse.ArgumentParser(description="Ensure unique sequences in alignments.")
    parser.add_argument('--fasta', type=str, required=True, help="Input FASTA file")
    parser.add_argument('--prefix', type=str, required=True, help="Output file prefix")
    args = parser.parse_args()

    ensure_unique_sequences(args.fasta, args.prefix)

if __name__ == "__main__":
    main()
