import argparse
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages
import heapq
import sys
import math

def parse_arguments():
    parser = argparse.ArgumentParser(description="Generate heatmaps from Selfish output")
    parser.add_argument("diff_file", help="Tab-delimited file containing differential bins")
    parser.add_argument("chrom_sizes", help="Tab-delimited chromosome sizes file")
    parser.add_argument("resolution", type=int, help="Binning resolution")
    parser.add_argument("--sample_name", help="Sample name (optional)")
    parser.add_argument("--gff_file", help="GFF file (optional)")
    parser.add_argument("--gene_numbers", help="Text file containing gene annotation numbers (optional)")
    return parser.parse_args()

def load_data(diff_file, chrom_sizes, resolution):
    diff_data = pd.read_csv(diff_file, sep='\t')
    chrom_sizes = pd.read_csv(chrom_sizes, sep='\t', header=None, names=['chr', 'length'])
    chrom_sizes['length'] = np.ceil(chrom_sizes['length'] / resolution).astype(int)
    return diff_data, chrom_sizes

def get_gene_locations(gff_file, gene_numbers_file, chrom_sizes, resolution):
    if gff_file is None:
        return None
    
    gff = pd.read_csv(gff_file, sep = '\t', header = None)

    if gene_numbers_file is None:
        genes = gff[gff[2] == 'gene']
    else:
        with open(gene_numbers_file, 'r') as f:
            gene_numbers = set(line.strip() for line in f)
            genes = gff.filter(lambda x: x[2] == 'gene' and x.attrs.get('ID', '').split(':')[-1] in gene_numbers)

    gene_locations = {}
    for _, row in chrom_sizes.iterrows():
        chrom, length = row['chr'], row['length']
        chrom_genes = genes[genes[0] == chrom]
        bins = [(int(min(gene[4] // resolution + 1, length)), 
                 int(min(gene[5] // resolution + 1, length))) 
                for gene in chrom_genes.itertuples()]
        unique_bins = list(set(bins))
        gene_locations[chrom] = unique_bins
    
    return gene_locations

def create_heatmap(diff_data, chrom, bins, resolution, pdf, gene_locations=None):
    matrix = np.zeros((bins, bins))
    chrom_data = diff_data[(diff_data['chr1'] == chrom) & (diff_data['chr2'] == chrom)]
    
    for _, row in chrom_data.iterrows():
        i, j = int(row['bin1'] / resolution), int(row['bin2'] / resolution)
        if i < bins and j < bins:
            matrix[i, j] = matrix[j, i] = row['log2FC']

    flat_matrix = np.abs(matrix.flatten())
    flat_matrix_filtered = flat_matrix[flat_matrix != 0]
    top10percent = int(np.ceil(len(flat_matrix_filtered) * 0.05))
    maxColor = min(heapq.nlargest(top10percent, flat_matrix_filtered))

    fig, ax = plt.subplots()
    plt.imshow(matrix, cmap='RdBu_r', aspect='equal', origin='lower', vmin=-maxColor, vmax=maxColor, extent=[0, bins, 0, bins])
    
    plt.title(f"Chromosome {chrom}", pad=12)
    ax.set_xlabel("Genomic Position (kb)")
    
    if gene_locations:
        chrom_genes = gene_locations.get(chrom, [])
        chrom_genes = sorted(set([start+0.5 for start, end in chrom_genes]))

        secaxx = ax.secondary_xaxis('top')
        secaxx.set_xticks(chrom_genes)
        secaxx.tick_params(width=1, length=8, color='#009e00', labeltop=False)
        secaxy = ax.secondary_yaxis('right')
        secaxy.set_yticks(chrom_genes)
        secaxy.tick_params(width=1, length=8, color='#009e00', labelright=False)

    plt.colorbar(label='log2FC')

    plt.tight_layout()
    pdf.savefig()
    plt.close()


def create_inter_heatmap(diff_data, chrom_sizes, resolution, pdf):
    df_inter = diff_data[diff_data.chr1 != diff_data.chr2]
    df_inter.loc[:, ['bin1', 'bin2']] = df_inter[['bin1', 'bin2']] // resolution + 1

    genome_length = chrom_sizes.length.sum()
    matrix = np.zeros((genome_length, genome_length))

    chrom_sizes['start'] = (chrom_sizes.length).cumsum()
    chrom_sizes['start'] = chrom_sizes['start'].shift(1, fill_value=0)

    df_inter_merged = pd.merge(df_inter, chrom_sizes, left_on='chr1', right_on='chr', how='left')
    df_inter_merged['bin1'] += df_inter_merged['start'].fillna(0)

    df_inter_merged = pd.merge(df_inter_merged, chrom_sizes, left_on='chr2', right_on='chr', suffixes=('_chr1', '_chr2'), how='left')
    df_inter_merged['bin2'] += df_inter_merged['start_chr2'].fillna(0)

    df_final = df_inter_merged.drop(columns=['chr_chr1', 'chr_chr2', 'start_chr1', 'start_chr2', 'length_chr1', 'length_chr2'])

    for row in df_final.itertuples():
        matrix[int(row.bin1) - 1][int(row.bin2) - 1] = row.log2FC
        matrix[int(row.bin2) - 1][int(row.bin1) - 1] = row.log2FC

    flat_matrix = np.abs(matrix.flatten())
    flat_matrix_filtered = flat_matrix[flat_matrix != 0]
    top10percent = int(np.ceil(len(flat_matrix_filtered) * 0.05))
    maxColor = min(heapq.nlargest(top10percent, flat_matrix_filtered))

    fig, ax = plt.subplots()
    plt.imshow(matrix, cmap='RdBu_r', aspect='equal', origin='lower', vmin=-maxColor, vmax=maxColor, extent=[0, genome_length, 0, genome_length])

    plt.title('Differential interchromosomal interactions', pad=12)
    
    ax.set_xticks(chrom_sizes['start'], minor=True)
    ax.set_yticks(chrom_sizes['start'], minor=True)
    
    ax.tick_params(which='minor', length=0, labeltop=False)
    ax.tick_params(which='minor', length=0, labeltop=False)
    
    ax.xaxis.grid(True, linestyle='dotted', linewidth=0.5, color='black', which='minor')
    ax.yaxis.grid(True, linestyle='dotted', linewidth=0.5, color='black', which='minor')

    plt.colorbar(label='log2FC')

    plt.tight_layout()
    pdf.savefig()
    plt.close()

def main():
    args = parse_arguments()
    
    diff_data, chrom_sizes = load_data(args.diff_file, args.chrom_sizes, args.resolution)
    gene_locations = get_gene_locations(args.gff_file, args.gene_numbers, chrom_sizes, args.resolution)

    pdf = PdfPages(f"{args.sample_name}_heatmaps.pdf" if args.sample_name else "heatmaps.pdf")

    for _, row in chrom_sizes.iterrows():
        chrom, bins = row['chr'], row['length']
        create_heatmap(diff_data, chrom, bins, args.resolution, pdf, gene_locations)

    create_inter_heatmap(diff_data, chrom_sizes, args.resolution, pdf)
    pdf.close()

if __name__ == "__main__":
    main()