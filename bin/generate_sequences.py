import pyvolve
import os
import argparse

# set a mutation rate
mutation_rates = {
    "AC": 0.039,
    "AG": 0.310,
    "AT": 0.123,
    "CA": 0.140,
    "CG": 0.022,
    "CT": 3.028,
    "GA": 0.747,
    "GC": 0.113,
    "GT": 2.953,
    "TA": 0.056,
    "TC": 0.261,
    "TG": 0.036
}

def generate_sequences(tree_file, out_dir, genome_length, mutation_rates):
    # Ensure output directory exists
    os.makedirs(out_dir, exist_ok=True)
    
    # Read the tree
    tree = pyvolve.read_tree(file=tree_file)
    
    # Define the model with gamma-distributed rate categories
    rate_matrix = {"AC": 0.039, "AG": 0.310, "AT": 0.123, "CA": 0.140, "CG": 0.022, "CT": 3.028,
                   "GA": 0.747, "GC": 0.113, "GT": 2.953, "TA": 0.056, "TC": 0.261, "TG": 0.036}
    model = pyvolve.Model("nucleotide", {"mu": rate_matrix})
    
    # Gamma-distributed rate categories
    category_rates = [1.0, 1.1, 1.5, 2.0]
    category_probs = [0.9, 0.01, 0.04, 0.05]
    
    # Normalize category rates to sum to 1
    total_rate = sum(category_rates)
    category_rates = [rate / total_rate for rate in category_rates]
    
    model = pyvolve.Model("nucleotide", {"mu": mutation_rates, "rate_categories": category_rates, "category_weights": category_probs})
    
    # Define the partition and evolver
    partition = pyvolve.Partition(models=model, size=genome_length)
    evolver = pyvolve.Evolver(tree=tree, partitions=partition)
    
    # Generate the sequences
    out_path = os.path.join(out_dir, os.path.basename(tree_file).replace(".tre", ""))
    evolver(ratefile=f"{out_path}_ratefile.txt", infofile=f"{out_path}_infofile.txt", seqfile=f"{out_path}_seqfile.fasta", write_anc=True)

def main():
    parser = argparse.ArgumentParser(description="Generate sequences from a phylogenetic tree with specific mutation rates.")
    parser.add_argument('--tree', type=str, required=True, help="Input tree file")
    parser.add_argument('--out', type=str, required=True, help="Output directory for sequences")
    parser.add_argument('--genome_length', type=int, required=True, help="Length of the genome sequences to generate")
    args = parser.parse_args()

    generate_sequences(args.tree, args.out, args.genome_length, mutation_rates)

if __name__ == "__main__":
    main()
