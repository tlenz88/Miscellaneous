from Bio import SeqIO

def filter_chromosomes(input_fasta, output_fasta, chromosome_prefixes):
    """
    Filters a FASTA file to keep only whole chromosomes.

    Parameters:
    - input_fasta: Path to the input FASTA file.
    - output_fasta: Path to the output FASTA file with only whole chromosomes.
    - chromosome_prefixes: A list of prefixes that identify whole chromosomes, e.g., ["chr", "chromosome"].
    """
    with open(input_fasta, "r") as infile, open(output_fasta, "w") as outfile:
        for record in SeqIO.parse(infile, "fasta"):
            # Check if the sequence ID starts with any of the chromosome prefixes
            if any(record.id.startswith(prefix) for prefix in chromosome_prefixes):
                SeqIO.write(record, outfile, "fasta")

if __name__ == "__main__":
    import sys
    input_fasta = sys.argv[1]
    output_fasta = sys.argv[2]
    
    # Adjust these prefixes based on how whole chromosomes are named in your FASTA file
    chromosome_prefixes = ["chr", "chromosome", "NC_"]  # Common prefixes for whole chromosomes

    filter_chromosomes(input_fasta, output_fasta, chromosome_prefixes)
