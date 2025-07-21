import sys
import argparse
import pysam
import pandas as pd
import numpy as np
from scipy import stats
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score
import matplotlib.pyplot as plt
import seaborn as sns
from statsmodels.nonparametric.smoothers_lowess import lowess
import re
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.gridspec as gridspec
import urllib.parse
from urllib.parse import unquote

def get_read_counts(bam_file):
    counts = {}
    with pysam.AlignmentFile(bam_file, "rb") as bam:
        for read in bam.fetch():
            chrom = read.reference_name
            pos = read.reference_start
            if chrom not in counts:
                counts[chrom] = {}
            counts[chrom][pos] = counts[chrom].get(pos, 0) + 1
    return counts

def normalize_coverage(chip_counts, input_counts):
    normalized = {}
    for chrom in chip_counts:
        normalized[chrom] = {}
        for pos in chip_counts[chrom]:
            chip_count = chip_counts[chrom].get(pos, 0)
            input_count = input_counts[chrom].get(pos, 0)
            if input_count > 0:
                normalized[chrom][pos] = np.log2(chip_count / input_count)
            else:
                normalized[chrom][pos] = 0
    return normalized

def process_gff(gff_file):
    gff_data = pd.read_csv(gff_file, sep='\t', comment='#', header=None,
                           names=['chromosome', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute'])
    
    gff_genes = gff_data[gff_data['feature'].isin(['gene', 'protein_coding_gene'])]
    
    def extract_attribute(attr, field):
        match = re.search(f'{field}=([^;]+)', attr)
        if match:
            return unquote(match.group(1))  # Decode URL-encoded characters
        return None
    
    gff_genes['gene_id'] = gff_genes['attribute'].apply(lambda x: extract_attribute(x, 'ID') or 
                                                                  extract_attribute(x, 'gene_id') or 
                                                                  extract_attribute(x, 'GeneID'))
    gff_genes['Name'] = gff_genes['attribute'].apply(lambda x: extract_attribute(x, 'Name'))
    gff_genes['Description'] = gff_genes['attribute'].apply(lambda x: extract_attribute(x, 'description') or 
                                                                      extract_attribute(x, 'Description'))
    gff_genes['ebi_biotype'] = gff_genes['attribute'].apply(lambda x: extract_attribute(x, 'ebi_biotype') or 
                                                                      extract_attribute(x, 'biotype'))

    # Decode URL-encoded characters in the description
    gff_genes['Description'] = gff_genes['Description'].apply(lambda x: urllib.parse.unquote(x) if pd.notna(x) else x)
    
    return gff_genes[['chromosome', 'start', 'end', 'strand', 'gene_id', 'Name', 'Description', 'ebi_biotype']]

def create_gene_families(deseq2_file, gff_data):
    deseq2_data = pd.read_csv(deseq2_file, sep='\t')
    
    merged_data = pd.merge(deseq2_data, gff_data[['gene_id', 'Name', 'Description']], on='gene_id', how='left')
    
    def assign_family(row):
        if pd.notna(row['Name']):
            return row['Name']
        elif pd.notna(row['Description']):
            # Split the description and take the first part (before the comma)
            return row['Description'].split(',')[0]
        else:
            return 'Unknown'
    
    merged_data['gene_family'] = merged_data.apply(assign_family, axis=1)
    
    gene_families = merged_data.groupby('gene_family')['gene_id'].apply(list).to_dict()
    
    return gene_families, merged_data[['gene_id', 'gene_family']]

def bin_region(counts, start, end, n_bins):
    region = [counts.get(pos, 0) for pos in range(start, end)]
    return np.array_split(region, n_bins)

def process_gene(gene, gff, counts):
    gene_info = gff[gff['gene_id'] == gene]
    chrom = gene_info['chromosome'].iloc[0]
    start = gene_info['start'].iloc[0]
    end = gene_info['end'].iloc[0]

    upstream_gene = gff[(gff['chromosome'] == chrom) & (gff['end'] < start)].iloc[-1] if not gff[(gff['chromosome'] == chrom) & (gff['end'] < start)].empty else None
    downstream_gene = gff[(gff['chromosome'] == chrom) & (gff['start'] > end)].iloc[0] if not gff[(gff['chromosome'] == chrom) & (gff['start'] > end)].empty else None

    upstream = bin_region(counts[chrom], upstream_gene['end'] if upstream_gene is not None else start - 1000, start, 5)
    gene_body = bin_region(counts[chrom], start, end, 5)
    downstream = bin_region(counts[chrom], end, downstream_gene['start'] if downstream_gene is not None else end + 1000, 5)
    upstream_gene_body = bin_region(counts[chrom], upstream_gene['start'] if upstream_gene is not None else start - 2000, upstream_gene['end'] if upstream_gene is not None else start - 1000, 2)
    downstream_gene_body = bin_region(counts[chrom], downstream_gene['start'] if downstream_gene is not None else end + 1000, downstream_gene['end'] if downstream_gene is not None else end + 2000, 2)

    result = np.concatenate([upstream_gene_body, upstream, gene_body, downstream, downstream_gene_body])
    return result


def find_optimal_k(data, max_k):
    silhouette_scores = []
    inertias = []
    for k in range(2, max_k+1):
        kmeans = KMeans(n_clusters=k, random_state=42)
        kmeans.fit(data)
        silhouette_scores.append(silhouette_score(data, kmeans.labels_))
        inertias.append(kmeans.inertia_)
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
    ax1.plot(range(2, max_k+1), silhouette_scores)
    ax1.set_xlabel('Number of clusters')
    ax1.set_ylabel('Silhouette Score')
    ax1.set_title('Silhouette Score vs. Number of Clusters')
    
    ax2.plot(range(2, max_k+1), inertias)
    ax2.set_xlabel('Number of clusters')
    ax2.set_ylabel('Inertia')
    ax2.set_title('Elbow Curve')
    
    plt.tight_layout()
    plt.savefig('optimal_k_plots.png')
    plt.close()

def process_deseq2_output(deseq2_file, genes):
    deseq2_data = pd.read_csv(deseq2_file, sep='\t')
    deseq2_data = deseq2_data[deseq2_data['gene_id'].isin(genes)]
    deseq2_data = deseq2_data.sort_values('log2FoldChange')
    return deseq2_data

def generate_combined_figure(wt_binned, treat_binned, ordered_genes, gene_families, clusters, output_file):
    with PdfPages(output_file) as pdf:
        fig = plt.figure(figsize=(20, 15))
        gs = gridspec.GridSpec(1, 5, width_ratios=[1, 1, 1, 0.2, 1])

        # Column 1: WT heatmap
        ax1 = plt.subplot(gs[0])
        sns.heatmap(wt_binned[ordered_genes], cmap='YlOrRd', xticklabels=False, yticklabels=False, ax=ax1)
        ax1.set_title("Wild Type H3K9me3/Input")

        # Column 2: Treatment heatmap
        ax2 = plt.subplot(gs[1])
        sns.heatmap(treat_binned[ordered_genes], cmap='YlOrRd', xticklabels=False, yticklabels=False, ax=ax2)
        ax2.set_title("Treatment H3K9me3/Input")

        # Column 3: Differential heatmap
        ax3 = plt.subplot(gs[2])
        sns.heatmap(treat_binned[ordered_genes] - wt_binned[ordered_genes], cmap='RdBu_r', xticklabels=False, yticklabels=False, ax=ax3, center=0)
        ax3.set_title("Log2(Treatment/Control)")

        # Column 4: Gene families
        ax4 = plt.subplot(gs[3])
        family_colors = sns.color_palette("husl", len(set(gene_families)))
        family_color_dict = dict(zip(set(gene_families), family_colors))
        family_colors_mapped = [family_color_dict[family] for family in gene_families]
        ax4.barh(range(len(gene_families)), [1]*len(gene_families), color=family_colors_mapped)
        ax4.set_ylim(ax4.get_ylim()[::-1])  # Invert y-axis to match heatmaps
        ax4.set_title("Gene Families")
        ax4.axis('off')

        # Column 5: LOESS plot
        ax5 = plt.subplot(gs[4])
        for cluster in range(7):
            cluster_data_wt = wt_binned[np.array(ordered_genes)[clusters == cluster]]
            cluster_data_treat = treat_binned[np.array(ordered_genes)[clusters == cluster]]
            x = np.arange(cluster_data_wt.shape[1])
            y_wt = np.mean(cluster_data_wt, axis=0)
            y_treat = np.mean(cluster_data_treat, axis=0)
            
            lowess_wt = lowess(y_wt, x, frac=0.6)
            lowess_treat = lowess(y_treat, x, frac=0.6)
            
            ax5.plot(lowess_wt[:, 0], lowess_wt[:, 1], label=f'WT Cluster {cluster}')
            ax5.plot(lowess_treat[:, 0], lowess_treat[:, 1], label=f'Treat Cluster {cluster}', linestyle='--')

        ax5.set_xlabel('Genomic Position')
        ax5.set_ylabel('Log2(H3K9me3/Input)')
        ax5.set_title('LOESS Regression of H3K9me3 Distribution')
        ax5.legend(bbox_to_anchor=(1.05, 1), loc='upper left')

        plt.tight_layout()
        pdf.savefig(fig)
        plt.close()

def main(args):

    # Step 1: Get read counts
    wt_h3k9me3_counts = get_read_counts(args.wt_h3k9me3)
    treat_h3k9me3_counts = get_read_counts(args.treat_h3k9me3)
    wt_input_counts = get_read_counts(args.wt_input)
    treat_input_counts = get_read_counts(args.treat_input)

    # Step 2: Normalize coverage
    wt_normalized = normalize_coverage(wt_h3k9me3_counts, wt_input_counts)
    treat_normalized = normalize_coverage(treat_h3k9me3_counts, treat_input_counts)

    # Step 3: Process GFF file
    gff_data = process_gff(args.gff)

    # Step 4-5: Create gene families and bin genes and neighboring regions
    gene_families, gene_family_df = create_gene_families(args.deseq2, gff_data)
    genes = gene_family_df['gene_id'].unique()

    wt_binned = np.array([process_gene(gene, gff_data, wt_normalized) for gene in genes])
    treat_binned = np.array([process_gene(gene, gff_data, treat_normalized) for gene in genes])

    # Step 6: Cluster genes
    find_optimal_k(wt_binned, 15)
    kmeans = KMeans(n_clusters=7, random_state=42)
    clusters = kmeans.fit_predict(wt_binned)

    # Step 7-8: Read differential expression data and order genes
    deseq2_data = process_deseq2_output(args.deseq2, genes)
    ordered_genes = []
    for cluster in range(7):
        cluster_genes = deseq2_data[deseq2_data['gene_id'].isin(genes[clusters == cluster])]
        ordered_genes.extend(cluster_genes['gene_id'].tolist())

    # Save gene families to a file
    gene_family_df.to_csv("gene_families.txt", sep='\t', index=False)

    # Generate combined figure
    generate_combined_figure(wt_binned, treat_binned, ordered_genes, gene_family_df['gene_family'], clusters, "combined_analysis.pdf")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Analyze H3K9me3 ChIP-seq data")
    parser.add_argument("--wt_h3k9me3", required=True, help="Wild type H3K9me3 BAM file")
    parser.add_argument("--treat_h3k9me3", required=True, help="Treatment H3K9me3 BAM file")
    parser.add_argument("--wt_input", required=True, help="Wild type Input BAM file")
    parser.add_argument("--treat_input", required=True, help="Treatment Input BAM file")
    parser.add_argument("--deseq2", required=True, help="DESeq2 output file")
    parser.add_argument("--gff", required=True, help="GFF or GTF file")
    args = parser.parse_args()

    main(args)
