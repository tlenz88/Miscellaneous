import pandas as pd
import argparse

# Argument parser to handle input files and output file
parser = argparse.ArgumentParser(description="Compare two files and output matching rows.")
parser.add_argument("-f1", "--file1", required=True, help="Path to the first tab-delimited text file")
parser.add_argument("-f2", "--file2", required=True, help="Path to the second tab-delimited text file")
parser.add_argument("-o", "--output", required=True, help="Path to the output file")

args = parser.parse_args()

# Read the input files
file1_df = pd.read_csv(args.file1, sep="\t", header=None)
file2_df = pd.read_csv(args.file2, sep="\t", header=None)

# Extract the column for comparison
file2_col9 = file2_df.iloc[:, 8]  # Column 9 of file2

# Find the matching rows from file1 based on the first column
matching_df = file1_df[file1_df.iloc[:, 0].isin(file2_col9)]

# Write the matching rows to the output file
matching_df.to_csv(args.output, sep="\t", index=False, header=False)
