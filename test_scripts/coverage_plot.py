import sys
import pysam
import pandas as pd
import matplotlib.pyplot as plt
from collections import defaultdict
import numpy as np
from scipy import stats
import argparse
import os
import re
from statsmodels.formula.api import glm
from statsmodels.genmod.families import NegativeBinomial
from statsmodels.stats.multitest import multipletests

def parse_arguments():
    parser = argparse.ArgumentParser(description="Analyze BAM files and plot read counts.")
    parser.add_argument("--control", nargs="+", required=True, help="Control/wild-type condition name followed by BAM files.")
    parser.add_argument("--conditions", nargs="+", action='append', help="Condition name followed by BAM files. Use multiple times for multiple conditions.")
    parser.add_argument("--inputs", nargs="+", help="BAM files for inputs (optional)")
    parser.add_argument("--num_bins", type=int, help="Number of bins for each region", default=5)
    parser.add_argument("--flank_length", type=int, help="Flanking region length", default=1000)
    parser.add_argument("--gff", help="GFF file", required=True)
    parser.add_argument("--gene_list", help="Text file containing list of genes")
    return parser.parse_args()

def parse_gff(gff_file, gene_list=None):
    genes = {}
    gene_ids = set()
    if gene_list:
        with open(gene_list, 'r') as f:
            gene_ids = set([line.strip() for line in f])

    with open(gff_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            if len(fields) != 9:
                continue
            
            feature_type = fields[2]
            attributes = dict(item.split('=') for item in fields[8].split(';') if '=' in item)
            
            if feature_type in ['protein_coding_gene', 'gene', 'ncRNA_gene']:
                gene_id = attributes.get('ID')
                if not gene_list or gene_id in gene_ids:
                    genes[gene_id] = {
                        'chr': fields[0],
                        'start': int(fields[3]),
                        'end': int(fields[4]),
                        'strand': fields[6],
                        'type': feature_type,
                        'exons': []
                    }
            
            elif feature_type == 'exon':
                exon_id = attributes.get('ID', '')
                
                gene_id_match = re.search(r'(PF3D7_\d+)', exon_id)
                if gene_id_match:
                    gene_id = gene_id_match.group(1)
                    if gene_id in genes:
                        genes[gene_id]['exons'].append((int(fields[3]), int(fields[4])))

    return genes

def sort_exons_by_strand(exons, strand):
    if strand == '+':
        return sorted(exons)
    else:
        return sorted(exons, reverse=True)

def bin_gene_regions(genes, num_bins, flanking_length):
    gene_bins = defaultdict(lambda: defaultdict(list))
    bin_sizes = defaultdict(lambda: defaultdict(list))
    for gene, info in genes.items():
        chr = info['chr']
        start = info['start']
        end = info['end']
        strand = info['strand']
        exons = sort_exons_by_strand(info['exons'], strand)
        if len(exons) < 2:
            continue

        CDS = (start, end)
        if strand == '+':
            intron = (exons[0][1], exons[1][0])
            five_prime = (start - flanking_length, start)
            exon1 = exons[0]
            exon2 = exons[1]
            three_prime = (end, end + flanking_length)
        else:
            intron = (exons[1][0], exons[0][1])
            five_prime = (end, end + flanking_length)
            exon1 = exons[1]
            exon2 = exons[0]
            three_prime = (start - flanking_length, start)

        regions = {
            'five_prime': five_prime,
            'exon1': exon1,
            'intron': intron,
            'exon2': exon2,
            'three_prime': three_prime
        }

        # Count reads in bins
        for region_name, (region_start, region_end) in regions.items():
            region_len = region_end - region_start
            bin_size_adjusted = region_len // num_bins
            gene_bins[gene]['chr'] = chr
            gene_bins[gene][region_name] = [(region_start + i * bin_size_adjusted, region_start + (i + 1) * bin_size_adjusted) for i in range(num_bins)]
            for i in range(num_bins):
                bin_sizes[gene][f'{region_name}_bin{i+1}'] = bin_size_adjusted

    return gene_bins, bin_sizes

def count_reads_in_bins(bam_file, gene_bins, num_bins):
    read_counts_binned = defaultdict(lambda: defaultdict(int))
    bam = pysam.AlignmentFile(bam_file, 'rb')
    for gene, info in gene_bins.items():
        chr = info['chr']
        for region_name, bin_coordinates in info.items():
            if region_name == 'chr':
                continue
            else:
                for i, (bin_start, bin_end) in enumerate(bin_coordinates):
                    read_counts_binned[gene][f"{region_name}_bin{i + 1}"] = 0
                    for read in bam.fetch(chr, bin_start, bin_end):
                        read_counts_binned[gene][f"{region_name}_bin{i + 1}"] += 1
    total_reads = bam.mapped

    bam.close()
    return read_counts_binned, total_reads

def normalize_counts(read_counts, total_reads, bin_sizes):
    length_normalized_counts = defaultdict(lambda: defaultdict(int))
    for gene, bins in read_counts.items():
        length_normalized_counts[gene] = {}
        for bin, count in bins.items():
            RPM_normalized_count = (count / total_reads) * 1e6
            length_normalized_counts[gene][bin] = RPM_normalized_count / bin_sizes[gene][bin]
    return length_normalized_counts

def input_normalize(condition_counts, input_counts):
    normalized_counts = {}
    for region in condition_counts:
        normalized_counts[region] = condition_counts[region] - input_counts[region]
    return normalized_counts

def merge_replicates(sample_data):
    merged_data = defaultdict(lambda: defaultdict(int))
    for rep in sample_data:
        for gene in sample_data[rep]:
            for bin in sample_data[rep][gene]:
                if bin not in merged_data[gene]:
                    merged_data[gene][bin] = []
                merged_data[gene][bin].append(sample_data[rep][gene][bin])

    for gene in merged_data:
        for bin in merged_data[gene]:
            merged_data[gene][bin] = np.mean(merged_data[gene][bin])

    return merged_data

def calculate_mean_counts_per_bin(counts):
    bin_counts = {i: [] for i in range(len(next(iter(counts.values()))))}
    for gene, bins in counts.items():
        for i, (bin, count) in enumerate(bins.items()):
            bin_counts[i].append(count)
    mean_counts = {i: np.mean(bin_counts[i]) for i in bin_counts}
    return mean_counts

def plot_results(counts_dict, p_values_dict, output_file):
    fig, ax = plt.subplots(figsize=(16, 8))
    
    # Flatten the counts for each condition
    flattened_counts = {condition: np.concatenate(list(counts.values())) 
                        for condition, counts in counts_dict.items()}

    x = np.arange(len(next(iter(flattened_counts.values()))))
    colors = plt.cm.rainbow(np.linspace(0, 1, len(flattened_counts)))
    
    for (condition, counts), color in zip(flattened_counts.items(), colors):
        if counts.size > 0:  # Only plot if there are counts
            ax.plot(x, counts, label=condition, color=color, marker='o', markersize=4, linewidth=1)
        else:
            print(f"Warning: No data to plot for condition {condition}")
    
    # Calculate region boundaries
    region_boundaries = [0]
    region_names = []
    for region, region_counts in next(iter(counts_dict.values())).items():
        region_boundaries.append(region_boundaries[-1] + len(region_counts))
        region_names.append(region)
    
    # Add vertical lines between regions
    for boundary in region_boundaries[1:-1]:
        ax.axvline(x=boundary - 0.5, color='gray', linestyle='--', alpha=0.5)
    
    # Set x-ticks and labels
    ax.set_xticks([(start + end) / 2 for start, end in zip(region_boundaries[:-1], region_boundaries[1:])])
    ax.set_xticklabels(region_names, rotation=45, ha='right')
    
    ax.set_ylabel('Normalized read count')
    ax.set_title('Read count per bin')
    ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    
    # Highlight significant bins
    for condition, p_values in p_values_dict.items():
        significant_bins = []
        current_bin = 0
        for region, p_value in p_values.items():
            if np.isscalar(p_value):
                if p_value < 0.05:
                    significant_bins.extend(range(current_bin, current_bin + len(counts_dict[condition][region])))
            else:
                significant_bins.extend([current_bin + i for i, p in enumerate(p_value) if p < 0.05])
            current_bin += len(counts_dict[condition][region])
        
        if significant_bins:
            ax.scatter(significant_bins, flattened_counts[condition][significant_bins], 
                       color='red', s=40, zorder=3, label=f'{condition} significant')
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()

def calculate_p_values(condition1, condition2):
    p_values = {}
    for region in condition1:
        _, p_values[region] = stats.ttest_ind(condition1[region], condition2[region], axis=0)
    return p_values

def save_p_values(p_values, genes, binned_regions, output_file):
    with open(output_file, 'w') as f:
        f.write("P-values for each region:\n")
        for region, p_value in p_values.items():
            if np.isscalar(p_value):
                f.write(f"{region}: {p_value}\n")
            else:
                f.write(f"{region}: {np.mean(p_value)}\n")
        
        f.write("\nP-values for each bin:\n")
        for region, bins in binned_regions.items():
            if np.isscalar(p_values[region]):
                for i in range(len(bins)):
                    f.write(f"{region} bin {i}: {p_values[region]}\n")
            else:
                for i, p_value in enumerate(p_values[region]):
                    f.write(f"{region} bin {i}: {p_value}\n")
        
        f.write("\nP-values for each region within each gene:\n")
        for gene in genes:
            for region in p_values:
                if region.startswith(gene):
                    if np.isscalar(p_values[region]):
                        f.write(f"{gene} - {region}: {p_values[region]}\n")
                    else:
                        f.write(f"{gene} - {region}: {np.mean(p_values[region])}\n")
        
        f.write("\nP-values for each bin within each gene:\n")
        for gene in genes:
            for region, bins in binned_regions.items():
                if region.startswith(gene):
                    if np.isscalar(p_values[region]):
                        for i in range(len(bins)):
                            f.write(f"{gene} - {region} bin {i}: {p_values[region]}\n")
                    else:
                        for i, p_value in enumerate(p_values[region]):
                            f.write(f"{gene} - {region} bin {i}: {p_value}\n")

def NB_test(df1_cond1, df2_cond1, df1_cond2, df2_cond2):
    # Function to perform negative binomial test for a single bin
    def nb_test_for_bin(data, bin_col):
        try:
            print(f"Fitting model for {bin_col}...")
            model = glm(f"{bin_col} ~ condition", data=data, family=NegativeBinomial(alpha=1.0))
            results = model.fit()
            print(f"Results for {bin_col}: {results.summary()}")
            condition_key = 'condition[T.condition2]'
            if condition_key in results.pvalues:
                return results.pvalues[condition_key], results.params[condition_key], results.bse[condition_key]
            else:
                print(f"Warning: '{condition_key}' not found in results for {bin_col}")
                return np.nan, np.nan, np.nan
        except Exception as e:
            print(f"Error in bin {bin_col}: {str(e)}")
            return np.nan, np.nan, np.nan

    # Combine the DataFrames
    df_cond1 = pd.concat([df1_cond1, df2_cond1]).reset_index(drop=True)
    df_cond2 = pd.concat([df1_cond2, df2_cond2]).reset_index(drop=True)

    df_cond1['condition'] = 'condition1'
    df_cond2['condition'] = 'condition2'

    combined_df = pd.concat([df_cond1, df_cond2]).reset_index(drop=True)

    # Convert condition to a categorical variable
    combined_df['condition'] = combined_df['condition'].astype('category')

    # Print combined DataFrame structure for debugging
    print("Combined DataFrame:")
    print(combined_df.head())
    print(combined_df.dtypes)

    # List of bin columns
    bin_columns = [col for col in combined_df.columns if col != 'condition']

    # Perform negative binomial test for each bin
    results = []
    for bin_col in bin_columns:
        p_value, coef, std_err = nb_test_for_bin(combined_df, bin_col)
        results.append((bin_col, p_value, coef, std_err))

    # Create results DataFrame
    results_df = pd.DataFrame(results, columns=['Bin', 'P-value', 'Coefficient', 'Std Error'])

    # Correct for multiple testing (ignoring NaN values)
    valid_p_values = results_df['P-value'].dropna()
    if not valid_p_values.empty:
        corrected_p_values = multipletests(valid_p_values, method='fdr_bh')[1]
        results_df.loc[valid_p_values.index, 'Adjusted P-value'] = corrected_p_values
    else:
        results_df['Adjusted P-value'] = np.nan

    # Sort results by adjusted p-value (NaN values at the end)
    results_df = results_df.sort_values('Adjusted P-value', na_position='last')

    print(results_df)

    # Identify significantly different bins (e.g., adjusted p-value < 0.05)
    significant_bins = results_df[results_df['Adjusted P-value'] < 0.05]
    print("\nSignificantly different bins:")
    print(significant_bins)

    # Calculate and print mean counts for significant bins in each condition
    if not significant_bins.empty:
        for bin_col in significant_bins['Bin']:
            mean_cond1 = df_cond1[bin_col].mean()
            mean_cond2 = df_cond2[bin_col].mean()
            print(f"\n{bin_col}:")
            print(f"Mean count in condition 1: {mean_cond1:.2f}")
            print(f"Mean count in condition 2: {mean_cond2:.2f}")

    # Print some diagnostic information
    print("\nDiagnostic Information:")
    print(f"Total bins: {len(bin_columns)}")
    print(f"Bins with valid p-values: {valid_p_values.count()}")
    print(f"Bins with errors or warnings: {len(bin_columns) - valid_p_values.count()}")

    # Check for zero variance in any bin
    zero_var_bins = combined_df[bin_columns].columns[combined_df[bin_columns].var() == 0].tolist()
    if zero_var_bins:
        print("\nWarning: The following bins have zero variance (all values are the same):")
        print(zero_var_bins)

    # Check for bins with all zeros
    all_zero_bins = combined_df[bin_columns].columns[combined_df[bin_columns].sum() == 0].tolist()
    if all_zero_bins:
        print("\nWarning: The following bins have all zero values:")
        print(all_zero_bins)

def main():
    args = parse_arguments()
    
    genes = parse_gff(args.gff, args.gene_list)
    gene_bins, bin_sizes = bin_gene_regions(genes, args.num_bins, args.flank_length)

    sample_counts = {}
    sample_totals = {}

    # Process control condition
    control_name = args.control[0]
    control_bams = args.control[1:]
    sample_counts[control_name] = {}
    sample_totals[control_name] = {}
    for i, bam_file in enumerate(control_bams):
        sample_counts[control_name][f'{control_name}_rep{i+1}'], sample_totals[control_name][f'{control_name}_rep{i+1}'] = count_reads_in_bins(bam_file, gene_bins, args.num_bins)

    # Process other conditions
    for condition in args.conditions:
        condition_name = condition[0]
        condition_bams = condition[1:]
        sample_counts[condition_name] = {}
        sample_totals[condition_name] = {}
        for i, bam_file in enumerate(condition_bams):
            sample_counts[condition_name][f'{condition_name}_rep{i+1}'], sample_totals[condition_name][f'{condition_name}_rep{i+1}'] = count_reads_in_bins(bam_file, gene_bins, args.num_bins)

    # Initialize lists to store DataFrames for each sample
    sample1_dfs = []
    sample2_dfs = []
    sample_list = list(sample_counts.values())
    # Iterate over the samples and replicates to create DataFrames
    for sample_index, samples in enumerate(sample_counts):
        for rep, genes in sample_counts[samples].items():
            # Convert the genes dictionary to a DataFrame
            df = pd.DataFrame.from_dict(genes, orient='index')
            if sample_index == 0:
                sample1_dfs.append(df)
            elif sample_index == 1:
                sample2_dfs.append(df)
    NB_test(sample1_dfs[0], sample1_dfs[1], sample2_dfs[0], sample2_dfs[1])
    sys.exit()

    normalized_sample_counts = {}
    for sample, reps in sample_counts.items():
        normalized_sample_counts[sample] = {}
        for rep, genes in reps.items():
            normalized_sample_counts[sample][rep] = normalize_counts(sample_counts[sample][rep], sample_totals[sample][rep], bin_sizes)

    merged_sample_counts = {}
    for sample, reps in normalized_sample_counts.items():
        if len(sample) > 1:
            merged_sample_counts[sample] = merge_replicates(normalized_sample_counts[sample])

    # Process input if available
    if args.inputs:
        input_counts = {}
        input_totals = {}
        for input in args.inputs:
            input_name = input[0]
            input_bams = input[1:]
            input_counts[input_name] = {}
            input_totals[input_name] = {}
            for i, bam_file in enumerate(input_bams):
                input_counts[input_name][f'{input_name}_rep{i+1}'], input_totals[input_name][f'{input_name}_rep{i+1}'] = count_reads_in_bins(bam_file, gene_bins, args.num_bins)

        normalized_input_counts = {}
        for sample, reps in input_counts.items():
            normalized_input_counts[sample] = {}
            for rep, genes in reps.items():
                normalized_input_counts[sample][rep] = normalize_counts(input_counts[sample][rep], input_totals[sample][rep], bin_sizes)

        merged_input_counts = {}
        for sample, reps in normalized_input_counts.items():
            if len(sample) > 1:
                merged_input_counts[sample] = merge_replicates(normalized_input_counts[sample])
        
        for sample in merged_input_counts:
            if sample in merged_sample_counts:
                merged_sample_counts[sample] = input_normalize(merged_sample_counts[sample], merged_input_counts[sample])

    sample_bin_counts = {}
    for sample in merged_sample_counts:
        sample_bin_counts[sample] = calculate_mean_counts_per_bin(merged_sample_counts[sample])

    # Calculate p-values and save results for each condition vs control
    p_values_dict = {}
    for sample in sample_bin_counts:
        if sample not in control_name:
            p_values = calculate_p_values(sample_bin_counts[control_name], sample_bin_counts[sample])
            p_values_dict[sample] = p_values
            save_p_values(p_values, genes, bin_sizes, f"p_values_{sample}_vs_{control_name}.txt")
    sys.exit()
    # Plot results
    plot_results(sample_bin_counts, p_values_dict, "results_plot.pdf")

if __name__ == "__main__":
    main()
