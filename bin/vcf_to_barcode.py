#!/usr/bin/env python
import argparse
import random
import pysam
import pandas as pd

def main():
    parser = argparse.ArgumentParser(description='Convert VCF to mutation table format.')
    parser.add_argument('vcf_file', help='Input VCF file')
    parser.add_argument('output_file', help='Output table file')
    parser.add_argument('--force_unique', action='store_true', help='Ensure no identical rows in the final table')
    args = parser.parse_args()

    # Load VCF file
    vcf = pysam.VariantFile(args.vcf_file)

    # Collect mutations and sample data
    mutations = []
    sample_names = vcf.header.samples
    data = {sample: [] for sample in sample_names}

    for record in vcf:
        mutation = f"{record.ref}{record.pos}{record.alts[0]}"
        mutations.append(mutation)
        for sample in sample_names:
            genotype = record.samples[sample]['GT']
            data[sample].append(genotype[0])
    # Create DataFrame
    df = pd.DataFrame(data, index=mutations).T

    if args.force_unique:
        ensure_unique_rows(df)

    # Save to output file
    df.to_csv(args.output_file, sep='\t')

def ensure_unique_rows(df):
    while df.duplicated(keep=False).any():
        duplicates = df[df.duplicated(keep=False)]
        for idx in duplicates.index:
            row = df.loc[idx]
            pos = random.randint(0, len(row) - 1)
            row.iloc[pos] = 1 - row.iloc[pos] # Flip random position
            df.loc[idx] = row

if __name__ == "__main__":
    main()
