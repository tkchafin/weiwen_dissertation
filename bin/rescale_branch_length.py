#!/usr/bin/env python
import os
from Bio import Phylo
import argparse

def scale_tree_branches_biopython(tree_file, output_file, scale_factor):
    # Load tree
    tree = Phylo.read(tree_file, "newick")

    # Rescale the branch length
    for clade in tree.find_clades():
        if clade.branch_length:
            clade.branch_length *= scale_factor

    # Write the scaled tree to output file
    Phylo.write(tree, output_file, "newick")

def main():
    # Setup command line argument parsing
    parser = argparse.ArgumentParser(description="Scale branch lengths of a phylogenetic tree.")
    parser.add_argument('--input_tree', type=str, required=True, help="Input Newick tree file")
    parser.add_argument('--scale', type=float, required=True, help="Scale factor for branch lengths")
    parser.add_argument('--output_dir', type=str, required=True, help="Directory where output Newick tree file will be saved")
    args = parser.parse_args()

    # Ensure output directory exists
    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)

    # Generate output filename with scale factor
    base_name = os.path.splitext(os.path.basename(args.input_tree))[0]  # Remove the extension
    scaled_filename = f"{base_name}_scale{args.scale}.tre"
    output_file = os.path.join(args.output_dir, scaled_filename)

    # Process the tree
    scale_tree_branches_biopython(args.input_tree, output_file, args.scale)

if __name__ == "__main__":
    main()
