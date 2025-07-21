import argparse
import pysam
import pandas as pd
import numpy as np
from scipy import stats
from collections import defaultdict

def parse_arguments():
    parser = argparse.ArgumentParser(description="Correlate ChIP-seq difference with differential expression.")
    parser.add_argument("--control", required=True, help="Control BAM file")
    parser.add_argument("--conditional", required=True, help="Conditional BAM file")
    parser.add_argument("--deseq2", required=True, help="DESeq2 output file")
    parser.add_argument("--gtf", required=True, help="GTF file for gene coordinates")
    parser.add_argument("--window", type=int, default=1000, help="Window size around TSS (default: 1000)")
    return parser.parse_args()

def read_gtf(gtf_file):
    genes = {}
    with open(gtf_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            if fields[2] == 'transcript':
                chrom = fields[0]
                start = int(fields[3])
                end = int(fields[4])
                strand = fields[6]
                attributes = dict(item.strip().split(' ') for item in fields[8].split(';') if item)
                gene_id = attributes.get('gene_id', '').strip('"')
                if gene_id:
                    genes[gene_id] = (chrom, start, end, strand)
    return genes

def count_reads(bam_file, genes, window):
    counts = defaultdict(int)
    with pysam.AlignmentFile(bam_file, "rb") as bam:
        for gene, (chrom, start, end, strand) in genes.items():
            if strand == '+':
                region_start = max(0, start - window)
                region_end = start + window
            else:
                region_start = max(0, end - window)
                region_end = end + window
            
            for read in bam.fetch(chrom, region_start, region_end):
                if not read.is_unmapped:
                    counts[gene] += 1
    return counts

def normalize_counts(counts, total_reads):
    return {gene: (count / total_reads) * 1e6 for gene, count in counts.items()}

def main():
    args = parse_arguments()

    # Read DESeq2 output
    deseq2_data = pd.read_csv(args.deseq2, sep='\t')
    
    # Read GTF file
    genes = read_gtf(args.gtf)

    # Count reads in control and conditional BAM files
    control_counts = count_reads(args.control, genes, args.window)
    conditional_counts = count_reads(args.conditional, genes, args.window)

    # Normalize counts
    total_control_reads = sum(control_counts.values())
    total_conditional_reads = sum(conditional_counts.values())
    norm_control_counts = normalize_counts(control_counts, total_control_reads)
    norm_conditional_counts = normalize_counts(conditional_counts, total_conditional_reads)

    # Calculate difference in ChIP-seq signal
    chip_diff = {gene: norm_conditional_counts.get(gene, 0) - norm_control_counts.get(gene, 0) 
                 for gene in set(norm_conditional_counts) | set(norm_control_counts)}

    # Prepare data for correlation
    corr_data = []
    for gene, log2fc in zip(deseq2_data['gene'], deseq2_data['log2FoldChange']):
        if gene in chip_diff:
            corr_data.append((log2fc, chip_diff[gene]))

    # Calculate correlation
    log2fc, chip_diff_values = zip(*corr_data)
    correlation, p_value = stats.pearsonr(log2fc, chip_diff_values)

    print(f"Pearson correlation: {correlation}")
    print(f"P-value: {p_value}")

    # Optional: create a scatter plot
    try:
        import matplotlib.pyplot as plt
        plt.figure(figsize=(10, 6))
        plt.scatter(log2fc, chip_diff_values, alpha=0.5)
        plt.xlabel('DESeq2 log2 Fold Change')
        plt.ylabel('ChIP-seq Signal Difference')
        plt.title(f'Correlation: {correlation:.3f} (p-value: {p_value:.3e})')
        plt.savefig('correlation_plot.png')
        print("Scatter plot saved as 'correlation_plot.png'")
    except ImportError:
        print("Matplotlib not installed. Skipping scatter plot generation.")

if __name__ == "__main__":
    main()
