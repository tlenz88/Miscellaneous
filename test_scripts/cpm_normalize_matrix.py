import pandas as pd
import numpy as np
import argparse

def cpm_normalize(input_file, output_file):
    # Read the .matrix file
    print(f"Reading input file: {input_file}")
    df = pd.read_csv(input_file, sep='\t', header=None, names=['bin1', 'bin2', 'count'])

    # Calculate total counts
    total_counts = df['count'].sum()
    print(f"Total counts: {total_counts}")

    # Perform CPM normalization
    print("Performing CPM normalization")
    df['cpm'] = (df['count'] / total_counts) * 1e6

    # Round CPM values to 6 decimal places
    df['cpm'] = df['cpm'].round(6)

    # Write the normalized matrix to a new file
    print(f"Writing normalized matrix to: {output_file}")
    df[['bin1', 'bin2', 'cpm']].to_csv(output_file, sep='\t', header=False, index=False)

    print("CPM normalization completed")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Perform CPM normalization on a HiC-Pro .matrix file")
    parser.add_argument("input_file", help="Path to the input .matrix file")
    parser.add_argument("output_file", help="Path to the output CPM-normalized .matrix file")
    args = parser.parse_args()

    cpm_normalize(args.input_file, args.output_file)