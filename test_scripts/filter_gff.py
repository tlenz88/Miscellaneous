import pandas as pd
import sys

# Function to filter the GFF file based on gene list
def filter_gff(gff_file, gene_list_file, output_file):
    # Load gene list
    with open(gene_list_file, 'r') as f:
        genes_to_keep = set(line.strip() for line in f if line.strip())

    # Open the GFF file and read line by line
    with open(gff_file, 'r') as infile, open(output_file, 'w') as outfile:
        for line in infile:
            if line.startswith("#"):  # Keep the header/comment lines
                continue
            else:
                fields = line.strip().split("\t")
                if len(fields) > 8:  # Ensure this is a proper GFF line
                    if fields[2] not in ["protein_coding_gene", "ncRNA_gene", "gene"]:
                        continue
                    attributes = fields[8]
                    # Extract gene information (assuming gene_id or similar is present)
                    gene_id = None
                    for attr in attributes.split(";"):
                        if attr.startswith("Name="):  # Check for 'gene=' in attributes
                            gene_id = attr.split("=")[-1]
                            break
                    
                    if gene_id and gene_id in genes_to_keep:
                        outfile.write(line)

# Example usage
gff_file = sys.argv[1]
gene_list_file = sys.argv[2]
output_file = sys.argv[3]

filter_gff(gff_file, gene_list_file, output_file)
