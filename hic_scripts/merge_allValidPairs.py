import os
import sys

def merge_allValidPairs(output_file, *input_files):
    """
    Merge multiple .allValidPairs files into a single output file.
    
    Parameters:
    - output_file: Path to the output merged file
    - input_files: Paths to input .allValidPairs files
    """
    with open(output_file, 'w') as outfile:
        for input_file in input_files:
            with open(input_file, 'r') as infile:
                for line in infile:
                    outfile.write(line)
    print(f"Merged {len(input_files)} .allValidPairs files into {output_file}.")


if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python merge_allValidPairs.py output_file input_file1 input_file2 ...")
        sys.exit(1)
    
    output_file = sys.argv[1]
    input_files = sys.argv[2:]
    
    # Check if input files exist
    for file in input_files:
        if not os.path.isfile(file):
            print(f"Error: {file} does not exist.")
            sys.exit(1)
    
    merge_allValidPairs(output_file, *input_files)
