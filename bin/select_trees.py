#!/usr/bin/env python
import os
import re
import argparse
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import dendropy
import numpy as np
from sklearn.linear_model import LinearRegression

# Function to parse tree file name
def parse_name(name):
    match = re.search(r's(\d+\.\d+)_r(\d+)', name)
    if match:
        return {'sd': float(match.group(1)), 'rep': int(match.group(2))}
    return {'sd': None, 'rep': None}

# Function to select and clean trees
def select_and_clean_trees(input_dir, output_dir, selection_file):
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    tree_files = [os.path.join(input_dir, f) for f in os.listdir(input_dir) if f.endswith('.tre')]
    trees = [(dendropy.Tree.get(path=f, schema="newick"), os.path.basename(f)) for f in tree_files]
    df = generate_tree_stats_df(trees)

    selected_tree_files = []
    if selection_file:
        try:
            with open(os.path.join(input_dir, selection_file), 'r') as file:
                selected_tree_files = [line.strip() for line in file if line.strip()]
        except FileNotFoundError:
            print(f"The file {selection_file} could not be read")
            return
    else:
        selected_tree_files = select_trees(df, out=output_dir)
    
    plot_tree_stats(df, selected_tree_files, out=output_dir)

    for filename_long in selected_tree_files:
        filename = os.path.basename(filename_long)
        file_path = os.path.join(input_dir, filename)
        if not os.path.exists(file_path):
            print(f"File {filename} listed in {selection_file} does not exist in the input directory.")
            continue

        with open(file_path, 'r') as infile:
            content = infile.read()
        cleaned_content = content.replace('[&R] ', '')

        with open(output_dir+"_"+filename, 'w') as outfile:
            outfile.write(cleaned_content)

# Function to generate tree stats dataframe
def generate_tree_stats_df(trees):
    data = []
    for tree, name in trees:
        colless = dendropy.calculate.treemeasure.colless_tree_imbalance(tree)
        sackin = dendropy.calculate.treemeasure.sackin_index(tree)
        treeness = dendropy.calculate.treemeasure.treeness(tree)
        gamma = dendropy.calculate.treemeasure.pybus_harvey_gamma(tree)
        b1 = dendropy.calculate.treemeasure.B1(tree)
        nbar = dendropy.calculate.treemeasure.N_bar(tree)
        parsed_params = parse_name(name)
        data.append({
            'SD': parsed_params['sd'],
            'Replicate': parsed_params['rep'],
            'Colless': colless,
            'Sackin': sackin,
            'Treeness': treeness,
            'Pybus&Harvey Gamma': gamma,
            'B1': b1,
            "Nbar": nbar,
            "Name": name
        })
    return pd.DataFrame(data)

def plot_tree_stats(df, selected_tree_files, out="out"):
    metrics = [col for col in df.columns if col not in ['SD', 'Replicate', 'Name', 'Predicted_Colless', 'Difference']]
    num_metrics = len(metrics)
    num_cols = 3
    num_rows = (num_metrics + num_cols - 1) // num_cols

    fig, axes = plt.subplots(nrows=num_rows, ncols=num_cols, figsize=(num_cols * 5, num_rows * 5))
    axes = axes.flatten()

    for i, metric in enumerate(metrics):
        sns.boxplot(x='SD', y=metric, data=df, ax=axes[i])
        sns.stripplot(x='SD', y=metric, data=df, color='k', ax=axes[i], alpha=0.6, jitter=True)
        axes[i].set_title(metric)
        axes[i].set_xlabel('Standard Deviation (SD)')
        axes[i].set_ylabel(metric)
        
        # Highlight selected trees
        selected_trees = df[df['Name'].isin(selected_tree_files)]
        for sd in selected_trees['SD'].unique():
            selected_subset = selected_trees[selected_trees['SD'] == sd]
            for _, row in selected_subset.iterrows():
                axes[i].scatter(x=[sd], y=[row[metric]], color='red', s=100, zorder=5, edgecolor='black', linewidth=1.5)

    for j in range(i + 1, len(axes)):
        fig.delaxes(axes[j])
    plt.tight_layout(rect=[0, 0, 1, 0.95])
    plt.savefig(f"{out}_stat_plots.pdf")
    plt.close()



def select_trees(df, out):
    # Sort sd and get the unique values
    sd_values = sorted(df['SD'].unique())
    selected_trees = []
    last_colless = None  # Store the colless value of the last selected tree

    for i, sd in enumerate(sd_values):
        # Select data for current 'sd' group
        current_group = df[df['SD'] == sd].sort_values(by='Colless', ascending=True).reset_index(drop=True)
        
        if i == 0:
            # For the first group: select the tree with min colless value
            selected_tree = current_group.iloc[0]
            last_colless = selected_tree['Colless']
        elif i == len(sd_values) - 1:
            # For the last group: select the tree with max colless value
            selected_tree = current_group.iloc[-1]
        else:
            # For other groups: find an appropriate tree
            median_index = len(current_group) // 2
            selected_tree = current_group.iloc[median_index]

            # Check if selected tree's colless value is greater than the last selected tree's colless value
            if selected_tree['Colless'] <= last_colless:
                # Find a valid tree that has a greater colless value than the last selected tree
                valid_trees = current_group.iloc[median_index:].query('Colless > @last_colless')
                if valid_trees.empty:
                    print(f"No valid tree found for sd={sd} that is greater than the last selected tree's colless value.")
                    continue  # Skip to the next group without updating last_colless
                selected_tree = valid_trees.iloc[0]

            last_colless = selected_tree['Colless']

        # Generate tree file name and add to the list
        tree_name = selected_tree['Name']
        selected_trees.append(tree_name)

    # Save names of selected trees
    n_selected = len(selected_trees)
    print(f"{n_selected} trees are selected.")
    with open(f"{out}_selected_trees.txt", 'w') as file:
        for tree in selected_trees:
            file.write(tree + '\n')

    return selected_trees



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Select and clean tree files.")
    parser.add_argument('--input', type=str, required=True, help="Input directory containing tree files")
    parser.add_argument('--out', type=str, required=True, help="Output directory for cleaned tree files")
    parser.add_argument('--select', type=str, default=None, help="File listing the trees to be processed")
    args = parser.parse_args()

    select_and_clean_trees(args.input, args.out, args.select)
