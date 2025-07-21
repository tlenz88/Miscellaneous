#!/usr/bin/env python3

import sys
import argparse
import os
import pysam
import matplotlib.pyplot as plt
from collections import defaultdict
import pandas as pd
from matplotlib.backends.backend_pdf import PdfPages
from scipy import stats
import re

def parse_args(args):
    parser = argparse.ArgumentParser()
    parser.add_argument('-c1',
                        '--condition1',
                        dest='condition1',
                        help='BAM files for conditional/experimental samples.',
                        nargs='+',
                        required=True)
    parser.add_argument('-i1',
                        '--input1',
                        dest='input1',
                        help='BAM files for conditional/experimental inputs.',
                        nargs='+',
                        required=False)
    parser.add_argument('-c2',
                        '--condition2',
                        dest='condition2',
                        help='BAM files for control/wild type samples.',
                        nargs='+',
                        required=True)
    parser.add_argument('-i2',
                        '--input2',
                        dest='input2',
                        help='BAM files for control/wild type inputs.',
                        nargs='+',
                        required=False)
    parser.add_argument('-o',
                        '--output',
                        dest='output',
                        help='Output PDF file name.',
                        required=False,
                        default='binned_gene_coverage.pdf')
    parser.add_argument('-b',
                        '--bins',
                        dest='bins',
                        help='Number of bins for each feature',
                        type=int,
                        default=5)
    parser.add_argument('-f',
                        '--flank',
                        help='Length of flanking regions',
                        type=int,
                        default=1000)
    parser.add_argument('-g',
                        '--gff',
                        dest='gff',
                        help='Gene data in GFF format. Coding regions are '
                             'used to determine which coordinates within BED '
                             'files to plot.',
                        required=True,
                        default=None)
    parser.add_argument('-l',
                        '--gene_list',
                        dest='gene_list',
                        help='List of genes to plot, separated by spaces, '
                             'commas, newlines or tabs. Genes will be '
                             'extracted from the provided GFF file.',
                        required=False,
                        default=None)
    parser.add_argument('-w',
                        '--whole_gene',
                        required=False,
                        action='store_true',
                        help='Plot the whole gene with flanking regions as a continuous plot',
                        default=False)
    parser.add_argument('-s',
                        '--summary',
                        dest='summary',
                        help='Output file for statistical summary.',
                        default='statistical_summary.txt')
    return parser.parse_args()

def parse_gff(gff_file, gene_list):
    gene_regions = {}
    if gene_list:
        gene_list = filter_gene_list(gene_list)
    with open(gff_file, 'r') as file:
        for line in file:
            if line.startswith('#'):
                continue
            cols = line.strip().split('\t')
            if cols[2] in ['gene', 'protein_coding_gene', 'ncRNA_gene']:
                gene_id = cols[8].split('ID=')[1].split(';')[0]
                if not gene_list or gene_id in gene_list:
                    chr = cols[0]
                    start = int(cols[3])
                    end = int(cols[4])
                    strand = cols[6]
                    gene_regions[gene_id] = {
                        'chr': chr,
                        'start': min(start, end),
                        'end': max(start, end),
                        'strand': strand,
                        'exons': []
                    }
            elif cols[2] == 'exon':
                exon_start = int(cols[3])
                exon_end = int(cols[4])
                try:
                    if cols[8].split('Parent=')[1].split('.')[1].split(';')[0] == '1':
                        parent_id = cols[8].split('Parent=')[1].split('.')[0]
                        if parent_id in gene_regions:
                            gene_regions[parent_id]['exons'].append((min(exon_start, exon_end), max(exon_start, exon_end)))
                except:
                    parent_id = cols[8].split('Parent=')[1].split('-')[0]
                    if parent_id in gene_regions:
                        gene_regions[parent_id]['exons'].append((min(exon_start, exon_end), max(exon_start, exon_end)))

    # Add exon count to each gene
    for gene, info in gene_regions.items():
        info['exon_count'] = len(info['exons'])

    return gene_regions

def filter_gene_list(gene_list):
    with open(gene_list, 'r') as gl:
        gene_list_str = gl.read()
    delimiters = "\t|,| |\n"
    gene_list = re.split(delimiters, gene_list_str)
    gene_list = [gene.strip() for gene in gene_list if gene.strip()]

    return gene_list

def get_sample_names(condition1, condition2):
    sample_list = []
    for sample in [condition1, condition2]:
        sample_name = os.path.basename(os.path.dirname(sample))
        if sample_name not in sample_list:
            sample_list.append(sample_name)
    return sample_list

def process_samples(args, condition_name, gene_regions, num_bins, flank):
    samples = getattr(args, condition_name, [])
    processed_samples = []
    
    for sample in samples:
        sample_name = os.path.splitext(os.path.basename(sample))[0]
        print(f"Counting and normalizing {sample_name} read counts")
        read_counts, bin_sizes, total_reads = count_reads_in_bins_bam(sample, gene_regions, num_bins, flank)
        norm_counts = normalize_per_million(read_counts, total_reads)
        normalized_counts = normalize_by_length(norm_counts, bin_sizes)
        processed_samples.append(normalized_counts)

    return processed_samples

def count_reads_in_bins_bam(bam_file, gene_regions, bin_count, flanking_length):
    read_counts_binned = defaultdict(lambda: defaultdict(int))
    bin_sizes = defaultdict(lambda: defaultdict(int))
    bam = pysam.AlignmentFile(bam_file, 'rb')

    for gene, info in gene_regions.items():
        chr = info['chr']
        start = min(info['start'], info['end'])
        end = max(info['start'], info['end'])
        strand = info['strand']
        exons = sorted(info['exons'], key=lambda x: x[0])

        if strand == '-':
            exons = exons[::-1]

        five_prime = (start - flanking_length, start) if strand == '+' else (end, end + flanking_length)
        three_prime = (end, end + flanking_length) if strand == '+' else (start - flanking_length, start)

        regions = {
            'five_prime': five_prime,
            'three_prime': three_prime
        }

        # Add exons and introns to regions
        for i, exon in enumerate(exons):
            exon_start, exon_end = min(exon), max(exon)
            regions[f'exon{i+1}'] = (exon_start, exon_end)
            if i < len(exons) - 1:
                if strand == '+':
                    intron_start = max(exon)
                    intron_end = min(exons[i+1])
                    regions[f'intron{i+1}'] = (min(intron_start, intron_end), max(intron_start, intron_end))
                elif strand == '-':
                    intron_start = min(exon)
                    intron_end = max(exons[i+1])
                    regions[f'intron{i+1}'] = (min(intron_start, intron_end), max(intron_start, intron_end))

        # Sort regions based on strand
        sorted_region_names = ['five_prime'] + \
                              [f'exon{i+1}' for i in range(len(exons))] + \
                              [f'intron{i+1}' for i in range(len(exons)-1)] + \
                              ['three_prime']

        if strand == '-':
            sorted_region_names = sorted_region_names[::-1]

        for region_name in sorted_region_names:
            for i in range(bin_count):
                bin_name = f"{region_name}_bin{i+1 if strand == '+' else bin_count-i}"
                read_counts_binned[gene][bin_name] = 0

        # Count reads in bins
        for region_name in sorted_region_names:
            region_start, region_end = regions[region_name]
            region_start, region_end = min(region_start, region_end), max(region_start, region_end)
            if region_start < 0:
                region_start = 0
            region_len = region_end - region_start
            bin_size_adjusted = max(region_len // bin_count, 1)
            bin_ranges = [(region_start + i * bin_size_adjusted, min(region_start + (i + 1) * bin_size_adjusted, region_end)) for i in range(bin_count)]

            if strand == '-':
                bin_ranges = bin_ranges[::-1]

            for i, (bin_start, bin_end) in enumerate(bin_ranges):
                bin_name = f"{region_name}_bin{i+1 if strand == '+' else bin_count-i}"
                bin_sizes[gene][bin_name] = bin_end - bin_start

                if bin_start >= bin_end:
                    continue
                try:
                    for read in bam.fetch(chr, bin_start, bin_end):
                        read_counts_binned[gene][bin_name] += 1
                except ValueError as e:
                    logging.warning(f"Warning: Error fetching reads for {gene}, {region_name}, bin {i+1}: {e}")
                    continue

    # Count total reads in the BAM file
    total_reads = bam.mapped
    bam.close()
    return read_counts_binned, bin_sizes, total_reads

def normalize_per_million(read_counts_binned, total_reads):
    """
    Normalize read counts by total reads and scale per-million (CPM).
    """
    normalized_counts = {}
    for gene, regions in read_counts_binned.items():
        normalized_counts[gene] = {}
        for region, count in regions.items():
            normalized_counts[gene][region] = (count / total_reads) * 1e6
    return normalized_counts

def normalize_by_length(normalized_counts, bin_sizes):
    """
    Normalize read counts by region lengths.
    """
    length_normalized_counts = {}
    for gene, regions in normalized_counts.items():
        length_normalized_counts[gene] = {}
        for region, count in regions.items():
            try:
                length_normalized_counts[gene][region] = count / bin_sizes[gene][region]
            except:
                length_normalized_counts[gene][region] = 0
    return length_normalized_counts

def remove_zero_count_genes(data):
    def is_non_zero(value):
        if isinstance(value, dict):
            return any(count != 0 for count in value.values())
        else:
            return value != 0

    return {
        sample: {
            gene: bin_counts
            for gene, bin_counts in genes.items()
            if is_non_zero(bin_counts)
        }
        for sample, genes in data.items()
    }

def merge_samples(*samples):
    """
    Merge the normalized counts by finding the mean across all provided samples.
    """
    if not samples:
        return defaultdict(lambda: defaultdict(float))

    merged_counts = defaultdict(lambda: defaultdict(float))
    sample_count = len(samples)

    for gene in samples[0]:
        for region in samples[0][gene]:
            total_count = 0
            valid_samples = 0
            for sample in samples:
                if gene in sample and region in sample[gene]:
                    total_count += sample[gene][region]
                    valid_samples += 1
            if valid_samples > 0:
                merged_counts[gene][region] = total_count / valid_samples

    return merged_counts

def input_subtraction(condition_counts, input_counts):
    input_subtracted_counts = defaultdict(lambda: defaultdict(float))
    for gene, regions in condition_counts.items():
        for region, count in regions.items():
            input_count = input_counts[gene][region]
            input_subtracted_counts[gene][region] = count - input_count
    return input_subtracted_counts

def read_counts_to_excel(all_read_counts_binned, output):
    combined_df = pd.DataFrame()

    for sample, genes in all_read_counts_binned.items():
        for gene, bins in genes.items():
            df = pd.DataFrame.from_dict(bins, orient='index', columns=[gene]).T
            combined_df = pd.concat([combined_df, df], axis=0)

    # Fill NaN values with 0 if you want to replace missing bin counts with 0
    combined_df = combined_df.fillna(0)

    # Save to Excel
    output_file = output.replace('.pdf', 'xlsx')
    combined_df.to_excel(output_file, sheet_name='GeneBinCounts')

def group_genes_by_exon_count(gene_regions):
    exon_groups = defaultdict(list)
    for gene, info in gene_regions.items():
        exon_count = info['exon_count']
        exon_groups[exon_count].append(gene)
    return exon_groups

def plot_read_counts(all_read_counts, samples, output, bin_count, gene_regions, plot_gene_body=False):
    plt.rcParams['font.family'] = 'serif'
    plt.rcParams['font.serif'] = ['Times New Roman'] + plt.rcParams['font.serif']

    exon_groups = group_genes_by_exon_count(gene_regions)

    def sort_regions(regions, strand, bin_count):
        sorted_names = ['five_prime']
        exon_count = int(sum(1 for name in regions if name.startswith('exon')) / bin_count)
        intron_count = exon_count - 1
        for i in range(1, exon_count + 1):
            sorted_names.append(f'exon{i}')
            if intron_count > 0 and i < exon_count:
                sorted_names.append(f'intron{i}')
        sorted_names.append('three_prime')
        return sorted_names

    pdf = PdfPages(output)

    for exon_count, genes in sorted(exon_groups.items()):
        plt.figure(figsize=(12, 8))
        
        for sample in all_read_counts.keys():
            counts = defaultdict(int)
            for gene in genes:
                if gene in all_read_counts[sample]:
                    strand = gene_regions[gene]['strand']
                    regions = sort_regions(all_read_counts[sample][gene].keys(), strand, bin_count)
                    for region in regions:
                        for i in range(bin_count):
                            bin_name = f"{region}_bin{i+1 if strand == '+' else bin_count-i}"
                            counts[bin_name] += all_read_counts[sample][gene][bin_name]

            x_labels = [bin for bin in counts.keys()]
            y_values = [count for count in counts.values()]

            plt.plot(range(1, len(x_labels) + 1), y_values, label=sample)

        plt.xlim(0, len(x_labels) + 0.5)
        ymin, ymax = plt.ylim()
        plt.ylim(0, ymax)

        ax = plt.gca()
        ax.spines['left'].set_position('zero')
        ax.spines['left'].set_color('black')
        ax.spines['left'].set_linewidth(1)
        ax.spines['bottom'].set_position('zero')
        ax.spines['bottom'].set_color('black')
        ax.spines['bottom'].set_linewidth(1)
        ax.spines['top'].set_color('none')
        ax.spines['right'].set_color('none')

        region_labels = ['5\' Region']
        if exon_count == 1:
            region_labels.append(f'Exon')
        elif exon_count > 1:
            for i in range(1, exon_count+1):
                region_labels.append(f'Exon {i}')
                if i < exon_count and exon_count > 2:
                    region_labels.append(f'Intron {i}')
                elif i < exon_count and exon_count == 2:
                    region_labels.append(f'Intron')
        region_labels.append('3\' Region')

        region_positions = [bin_count * (i + 1) + 0.5 for i in range(len(regions) - 1)]
        for pos in region_positions:
            plt.axvline(pos, color='grey', linestyle='--')

        plt.xticks([bin_count * (i + 0.5) for i in range(len(regions))], region_labels, rotation=45, ha='right')

        plt.ylabel('Normalized read count (CPM)')
        plt.title(f'Read counts per bin in gene regions (Genes with {exon_count} exons)')
        plt.legend()
        plt.tight_layout()
        pdf.savefig()
        plt.close()

    pdf.close()

def perform_statistical_tests(all_read_counts):
    samples = list(all_read_counts.keys())
    if len(samples) != 2:
        raise ValueError("This function expects exactly two samples for comparison")

    sample1, sample2 = samples
    condition1 = all_read_counts[sample1]
    condition2 = all_read_counts[sample2]

    # Collect all genes and bins from just these two conditions
    all_genes = set(condition1.keys()).union(set(condition2.keys()))
    all_bins = set()
    for gene_data in condition1.values():
        all_bins.update(gene_data.keys())
    for gene_data in condition2.values():
        all_bins.update(gene_data.keys())

    aggregated_counts1 = {bin: [] for bin in all_bins}
    aggregated_counts2 = {bin: [] for bin in all_bins}

    # Aggregate counts for each bin across all genes
    for gene in all_genes:
        bins1 = condition1.get(gene, {})
        bins2 = condition2.get(gene, {})

        for bin in all_bins:
            count1 = bins1.get(bin, 0)
            count2 = bins2.get(bin, 0)
            aggregated_counts1[bin].append(count1)
            aggregated_counts2[bin].append(count2)

    results = []

    # Perform rank-sum test for each bin
    for bin in all_bins:
        counts1 = aggregated_counts1[bin]
        counts2 = aggregated_counts2[bin]

        # Perform rank-sum test
        stat, p_value = stats.ranksums(counts1, counts2, alternative='two-sided')

        results.append({
            'bin': bin,
            'p-value': p_value,
            'statistic': stat
        })

    results_df = pd.DataFrame(results)
    return results_df

def write_statistical_summary(all_read_counts, gene_regions, output_file):
    """
    Generate a text file with summary statistical data.
    
    :param all_read_counts: Dictionary containing read counts for all samples
    :param gene_regions: Dictionary containing gene information
    :param output_file: Name of the output text file
    """
    # Group genes by exon count
    exon_groups = group_genes_by_exon_count(gene_regions)
    
    with open(output_file, 'w') as f:
        f.write("Statistical Summary\n")
        f.write("===================\n\n")
        
        for exon_count, genes in sorted(exon_groups.items()):
            f.write(f"Genes with {exon_count} exons:\n")
            f.write("-" * 30 + "\n")
            
            # Prepare data for this exon group
            group_data = {
                sample: {gene: counts for gene, counts in sample_data.items() if gene in genes}
                for sample, sample_data in all_read_counts.items()
            }
            
            # Perform statistical test
            results = perform_statistical_tests(group_data)
            
            # Write results
            for _, row in results.iterrows():
                bin_name = row['bin']
                p_value = row['p-value']
                statistic = row['statistic']
                f.write(f"Bin: {bin_name}\n")
                f.write(f"  p-value: {p_value:.6f}\n")
                f.write(f"  statistic: {statistic:.6f}\n")
            
            f.write("\n")

    print(f"Statistical summary has been written to {output_file}")

def main():
    args = parse_args(sys.argv[1:])

    sample_list = get_sample_names(args.condition1[0], args.condition2[0])

    gene_regions = parse_gff(args.gff, args.gene_list)
    num_bins = args.bins
    flank = args.flank

    condition_counts = {
        sample_list[0]: process_samples(args, 'condition1', gene_regions, num_bins, flank),
        sample_list[1]: process_samples(args, 'condition2', gene_regions, num_bins, flank)
    }

    merged_conditions = []
    for condition_name, samples in condition_counts.items():
        print(f"Merging replicates for {condition_name}")
        merged_conditions.append(merge_samples(*samples))

    if args.input1 and args.input2:
        input_counts = {
            'input1': process_samples(args, 'input1', gene_regions, num_bins, flank),
            'input2': process_samples(args, 'input2', gene_regions, num_bins, flank)
        }

        merged_inputs = []
        for input_name, samples in input_counts.items():
            print(f"Merging replicates for {input_name}")
            merged_inputs.append(merge_samples(*samples))

        all_read_counts_binned = {
            sample: input_subtraction(condition, input_)
            for sample, condition, input_ in zip(sample_list, merged_conditions, merged_inputs)
        }

    else:
        all_read_counts_binned = {
            sample: condition
            for sample, condition in zip(sample_list, merged_conditions)
        }

    plot_read_counts(all_read_counts_binned, sample_list, args.output, args.bins, gene_regions, args.whole_gene)

    summary_output = os.path.splitext(args.output)[0] + '_summary.txt'
    write_statistical_summary(all_read_counts_binned, gene_regions, summary_output)
    #results = perform_statistical_tests(all_read_counts_binned)
    #significant_bins = results[results['p-value'] < 0.05]
    #print(significant_bins)

if __name__ == "__main__":
    main()
