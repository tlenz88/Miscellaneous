import argparse
import csv

def parse_gff(gff_file, gene_list, upstream=0, downstream=0):
    regions = {}
    with open(gff_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            if len(fields) != 9:
                continue
            if fields[2] == 'gene':
                attributes = dict(item.split('=') for item in fields[8].split(';'))
                gene_id = attributes.get('ID', '').split(':')[-1]
                if gene_id in gene_list:
                    chrom = fields[0]
                    start = int(fields[3]) - upstream
                    end = int(fields[4]) + downstream
                    strand = fields[6]
                    regions[gene_id] = (chrom, start, end, strand)
    return regions

def write_bed(regions, output_file):
    with open(output_file, 'w') as f:
        for gene_id, (chrom, start, end, strand) in regions.items():
            f.write(f"{chrom}\t{start}\t{end}\t{gene_id}\t0\t{strand}\n")

def main():
    parser = argparse.ArgumentParser(description="Extract gene regions from GFF and create BED file")
    parser.add_argument("gff_file", help="Input GFF file")
    parser.add_argument("gene_list", help="File containing list of gene names, one per line")
    parser.add_argument("output_bed", help="Output BED file")
    parser.add_argument("--upstream", type=int, default=0, help="Number of base pairs to include upstream of each gene")
    parser.add_argument("--downstream", type=int, default=0, help="Number of base pairs to include downstream of each gene")
    args = parser.parse_args()

    # Read the list of genes
    with open(args.gene_list, 'r') as f:
        gene_list = set(line.strip() for line in f)

    # Parse the GFF file and extract regions
    regions = parse_gff(args.gff_file, gene_list, args.upstream, args.downstream)

    # Write the BED file
    write_bed(regions, args.output_bed)

    print(f"Extracted {len(regions)} regions to {args.output_bed}")

if __name__ == "__main__":
    main()
