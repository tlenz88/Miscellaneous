import sys
import pandas as pd
import numpy as np
import cooler
import scipy.stats as stats
from sklearn.preprocessing import StandardScaler
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import h5py
from scipy.sparse import csr_matrix
import gc

def load_and_normalize_expression(counts_file, metadata_file):
    """
    Load and normalize expression data using DESeq2-style normalization
    
    Parameters:
    counts_file: str, path to tab-delimited file with raw counts (genes x samples)
    metadata_file: str, path to sample metadata file with columns: sample_id, condition
    
    Returns:
    DataFrame with normalized mean expression per condition
    """
    # Load count data
    counts_df = pd.read_csv(counts_file, sep='\t', index_col=0)
    
    # Load metadata
    metadata = pd.read_csv(metadata_file, sep='\t')
    
    # Verify sample names match between counts and metadata
    assert all(counts_df.columns.isin(metadata['sample'])), "Sample names don't match between counts and metadata"
    
    # Calculate size factors (DESeq2-style)
    def geometric_mean(x):
        return np.exp(np.mean(np.log(x[x > 0])))
    
    # Calculate geometric means per gene
    geometric_means = counts_df.apply(geometric_mean, axis=1)
    
    # Calculate size factors
    size_factors = {}
    for sample in counts_df.columns:
        ratios = counts_df[sample] / geometric_means
        # Use median of ratios as size factor
        size_factors[sample] = np.median(ratios[~np.isnan(ratios) & ~np.isinf(ratios)])
    
    # Normalize counts
    normalized_counts = counts_df.copy()
    for sample in counts_df.columns:
        normalized_counts[sample] = counts_df[sample] / size_factors[sample]
    
    # Calculate mean normalized counts per condition
    condition_means = pd.DataFrame(index=normalized_counts.index)
    for condition in metadata['group'].unique():
        condition_samples = metadata[metadata['group'] == condition]['sample']
        condition_means[condition] = normalized_counts[condition_samples].mean(axis=1)
    
    return condition_means

def calculate_bin_interactions(cool_file, chunk_size=10000):
    """
    Calculate total interactions per bin using chunked processing
    
    Parameters:
    cool_file: str, path to .cool file
    chunk_size: int, number of rows to process at once
    
    Returns:
    DataFrame with bin information and total interactions
    """
    print(f"Processing Hi-C file: {cool_file}")
    c = cooler.Cooler(cool_file)
    
    # Get total number of bins
    n_bins = c.bins()[:].shape[0]
    
    # Initialize array for total interactions
    total_interactions = np.zeros(n_bins)
    
    # Process matrix in chunks
    for i in range(0, n_bins, chunk_size):
        end_idx = min(i + chunk_size, n_bins)
        print(f"Processing bins {i}-{end_idx} of {n_bins}")
        
        # Get chunk of matrix
        chunk = c.matrix(balance=True)[i:end_idx, :]
        
        # If chunk is sparse, convert to csr_matrix
        if isinstance(chunk, np.ndarray):
            chunk = csr_matrix(chunk)
        
        # Sum interactions for rows in chunk
        total_interactions[i:end_idx] = chunk.sum(axis=1).A1
        
        # Add column sums for this chunk (avoiding double-counting diagonal)
        col_sums = chunk.sum(axis=0).A1
        total_interactions += col_sums
        
        # Correct for double-counting diagonal elements
        diag_elements = chunk.diagonal()
        total_interactions[i:end_idx] -= diag_elements
        
        # Clear memory
        del chunk
        gc.collect()
    
    # Create DataFrame with bin information
    bins_df = c.bins()[:]
    bins_df['total_interactions'] = total_interactions
    
    return bins_df

def map_genes_to_bins(expression_df, hic_bins, genome_annotation):
    """
    Map genes to Hi-C bins using genomic coordinates
    """
    # Load gene annotations in chunks if file is large
    chunk_size = 100000
    gene_chunks = pd.read_csv(genome_annotation, sep='\t', chunksize=chunk_size)
    
    gene_bin_maps = []
    for genes_chunk in gene_chunks:
        gene_bin_map = pd.DataFrame()
        gene_bin_map['Gene_ID'] = genes_chunk['Gene_ID']
        gene_bin_map['bin_idx'] = ((genes_chunk['start'] + genes_chunk['end']) // 2) // hic_bins['resolution']
        gene_bin_maps.append(gene_bin_map)
    
    gene_bin_map = pd.concat(gene_bin_maps)
    
    # Merge with expression data
    mapped_data = gene_bin_map.merge(
        expression_df.reset_index().rename(columns={'index': 'Gene_ID'}),
        on='Gene_ID'
    )
    
    return mapped_data

def calculate_correlations(mapped_data, hic_bins, condition):
    """
    Calculate correlation between gene expression and interaction counts
    """
    # Process in chunks if dataset is large
    chunk_size = 10000
    pearson_corrs = []
    pearson_ps = []
    spearman_corrs = []
    spearman_ps = []
    
    for i in range(0, len(mapped_data), chunk_size):
        chunk = mapped_data.iloc[i:i + chunk_size]
        chunk_bins = hic_bins.loc[chunk['bin_idx']]
        
        # Prepare data
        scaler = StandardScaler()
        expr_scaled = scaler.fit_transform(chunk[[condition]].values)
        inter_scaled = scaler.fit_transform(chunk_bins[['total_interactions']].values)
        
        # Calculate correlations for chunk
        pearson = stats.pearsonr(expr_scaled.flatten(), inter_scaled.flatten())
        spearman = stats.spearmanr(expr_scaled.flatten(), inter_scaled.flatten())
        
        pearson_corrs.append(pearson[0])
        pearson_ps.append(pearson[1])
        spearman_corrs.append(spearman[0])
        spearman_ps.append(spearman[1])
    
    # Combine results (using Fisher's z-transformation for correlation coefficients)
    def fisher_combine(corrs, ps):
        z_scores = np.arctanh(corrs)  # Fisher's z-transformation
        mean_z = np.mean(z_scores)
        combined_corr = np.tanh(mean_z)  # Inverse transformation
        combined_p = stats.combine_pvalues(ps)[1]  # Fisher's method for p-values
        return combined_corr, combined_p
    
    final_pearson = fisher_combine(pearson_corrs, pearson_ps)
    final_spearman = fisher_combine(spearman_corrs, spearman_ps)
    
    return {
        'pearson': final_pearson,
        'spearman': final_spearman
    }

def plot_correlation(mapped_data, hic_bins, condition, output_prefix):
    """
    Create scatter plot of expression vs interactions using random sampling
    """
    # Sample data points for plotting (to avoid memory issues with large datasets)
    n_samples = min(10000, len(mapped_data))
    sampled_indices = np.random.choice(len(mapped_data), n_samples, replace=False)
    
    plot_data = mapped_data.iloc[sampled_indices].merge(
        hic_bins,
        left_on='bin_idx',
        right_index=True
    )
    
    plt.figure(figsize=(10, 8))
    sns.scatterplot(data=plot_data,
                   x=condition,
                   y='total_interactions',
                   alpha=0.5)
    
    plt.xlabel(f'Normalized Gene Expression ({condition})')
    plt.ylabel('Hi-C Interactions')
    plt.title(f'Gene Expression vs Hi-C Interactions\n{condition}')
    plt.savefig(f'{output_prefix}_{condition}_correlation_plot.png')
    plt.close()

def plot_condition_comparison(conditions_data, output_prefix):
    """
    Create comparison plots between conditions using random sampling
    """
    n_conditions = len(conditions_data)
    if n_conditions > 1:
        fig, axes = plt.subplots(1, n_conditions - 1, figsize=(6 * (n_conditions - 1), 5))
        if n_conditions == 2:
            axes = [axes]
        
        conditions = list(conditions_data.keys())
        n_samples = 10000  # Number of points to plot
        
        for i, ax in enumerate(axes):
            condition1 = conditions[i]
            condition2 = conditions[i + 1]
            
            # Sample points for plotting
            n_bins = len(conditions_data[condition1])
            sample_idx = np.random.choice(n_bins, n_samples, replace=False)
            
            ax.scatter(
                conditions_data[condition1]['total_interactions'].iloc[sample_idx],
                conditions_data[condition2]['total_interactions'].iloc[sample_idx],
                alpha=0.5
            )
            ax.set_xlabel(f'{condition1} Interactions')
            ax.set_ylabel(f'{condition2} Interactions')
            ax.set_title(f'Hi-C Interactions\n{condition1} vs {condition2}')
        
        plt.tight_layout()
        plt.savefig(f'{output_prefix}_condition_comparison.png')
        plt.close()

def main(counts_file, metadata_file, cool_files, genome_annotation, resolution, output_prefix):
    """
    Main function to run the correlation analysis
    """
    # Load and normalize expression data
    print("Loading and normalizing expression data...")
    expression_df = load_and_normalize_expression(counts_file, metadata_file)
    
    # Verify conditions match
    conditions = set(expression_df.columns)
    hic_conditions = set(cool_files.keys())
    assert conditions == hic_conditions, f"Conditions don't match between expression data ({conditions}) and Hi-C files ({hic_conditions})"
    
    # Process each condition
    results = {}
    conditions_data = {}
    
    for condition in conditions:
        print(f"Processing condition: {condition}")
        
        # Calculate bin interactions
        hic_bins = calculate_bin_interactions(cool_files[condition])
        conditions_data[condition] = hic_bins
        
        # Map genes to bins
        mapped_data = map_genes_to_bins(expression_df, hic_bins, genome_annotation)
        
        # Calculate correlations
        correlations = calculate_correlations(mapped_data, hic_bins, condition)
        results[condition] = correlations
        
        # Create correlation plot
        plot_correlation(mapped_data, hic_bins, condition, output_prefix)
    
    # Create comparison plots
    plot_condition_comparison(conditions_data, output_prefix)
    
    # Save results
    print("Saving results...")
    with open(f'{output_prefix}_correlation_results.txt', 'w') as f:
        for condition, corr in results.items():
            f.write(f'\nResults for {condition}:\n')
            f.write(f'Pearson correlation: {corr["pearson"][0]:.3f} (p={corr["pearson"][1]:.2e})\n')
            f.write(f'Spearman correlation: {corr["spearman"][0]:.3f} (p={corr["spearman"][1]:.2e})\n')
    
    print("Analysis complete!")

# Example usage
if __name__ == "__main__":
    counts_file = sys.argv[1]
    metadata_file = sys.argv[2]
    
    # Map conditions to their respective cool files
    cool_files = {
        'HFF': sys.argv[3],
        'ME49RFP_6hpi': sys.argv[4],
        'ME49RFP_24hpi': sys.argv[5]
    }
    
    genome_annotation = sys.argv[6]
    resolution = sys.argv[7]
    output_prefix = 'correlation_analysis'
    
    main(counts_file, metadata_file, cool_files, genome_annotation, resolution, output_prefix)