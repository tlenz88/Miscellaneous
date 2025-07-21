import pysam
import argparse
import os

def load_chromosome_mapping(mapping_file):
    """Load the chromosome name mapping from a tab-delimited text file."""
    mapping = {}
    with open(mapping_file, 'r') as f:
        for line in f:
            old_name, new_name = line.strip().split('\t')
            mapping[old_name] = new_name
    return mapping

def adjust_chromosome_names(input_bam, output_bam, mapping_file):
    """Adjust chromosome names in a BAM file based on a given mapping."""
    # Load chromosome mapping
    chrom_mapping = load_chromosome_mapping(mapping_file)
    
    # Open the input BAM file
    bam_in = pysam.AlignmentFile(input_bam, "rb")
    
    # Modify header to reflect the new chromosome names
    header = bam_in.header.to_dict()
    for seq in header['SQ']:
        if seq['SN'] in chrom_mapping:
            seq['SN'] = chrom_mapping[seq['SN']]
    
    # Open the output BAM file with modified header
    bam_out = pysam.AlignmentFile(output_bam, "wb", header=header)
    
    # Process each read in the input BAM file
    for read in bam_in.fetch(until_eof=True):
        # Adjust the reference name if it matches one in the mapping
        if bam_in.get_reference_name(read.reference_id) in chrom_mapping:
            new_name = chrom_mapping[bam_in.get_reference_name(read.reference_id)]
            read.reference_id = bam_in.get_tid(new_name)
        
        # Adjust the mate reference name if it matches one in the mapping
        if read.next_reference_id != -1 and bam_in.get_reference_name(read.next_reference_id) in chrom_mapping:
            new_name = chrom_mapping[bam_in.get_reference_name(read.next_reference_id)]
            read.next_reference_id = bam_in.get_tid(new_name)
        
        bam_out.write(read)

    # Close the files
    bam_in.close()
    bam_out.close()

    # Index the output BAM file
    pysam.index(output_bam)
    print(f"Chromosome names adjusted and saved to {output_bam}")

def main():
    parser = argparse.ArgumentParser(description="Adjust chromosome names in a BAM file.")
    parser.add_argument("-b", "--bam", required=True, help="Input BAM file")
    parser.add_argument("-m", "--mapping", required=True, help="Mapping file (tab-delimited) for chromosome name conversion")
    parser.add_argument("-o", "--output", help="Output BAM file (default: <input_basename>_updated.bam)")
    
    args = parser.parse_args()

    # If no output file is provided, use the default naming convention
    input_bam = args.bam
    mapping_file = args.mapping
    output_bam = args.output if args.output else os.path.splitext(input_bam)[0] + "_updated.bam"

    adjust_chromosome_names(input_bam, output_bam, mapping_file)

if __name__ == "__main__":
    main()
