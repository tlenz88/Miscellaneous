import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from sklearn.mixture import GaussianMixture
from scipy.stats import norm, percentileofscore
from scipy.optimize import fsolve
from sklearn.preprocessing import StandardScaler
from sklearn.cluster import KMeans
import sys
import argparse
import pysam
import re

def read_gff(gff_file):
    df = pd.read_csv(gff_file, sep="\t", header=None, comment="#")
    df = df[(df[2] == "protein_coding_gene") | (df[2] == "ncRNA_gene")]
    df.reset_index(drop=True, inplace=True)
    
    def extract_gene_id(attr):
        match = re.search(r"(?:ID|Name|gene_id)=?([^;\"\n]+)", attr)
        return match.group(1) if match else "N/A"
    
    bed = pd.DataFrame({
        "chrom": df[0],
        "start": df[3] - 1,
        "end": df[4],
        "name": df[8].apply(extract_gene_id),
        "score": 0,
        "strand": df[6]
    })
    return bed

def get_window(row, upstream, downstream):
    if upstream is None or downstream is None:
        return pd.Series([row["start"], row["end"]])
    
    strand = row.get("strand", "+")
    if strand == "+":
        tss = row["start"]
        win_start = max(0, tss - upstream)
        win_end = tss + downstream
    elif strand == "-":
        tss = row["end"]
        win_start = max(0, tss - downstream)
        win_end = tss + upstream
    else:
        tss = row["start"]
        win_start = max(0, tss - upstream)
        win_end = tss + downstream
    return pd.Series([win_start, win_end])

def count_reads(bam_file, bed_df, upstream, downstream):
    bam = pysam.AlignmentFile(bam_file, "rb")
    windows = bed_df.apply(get_window, axis=1, upstream=upstream, downstream=downstream)
    windows.columns = ["win_start", "win_end"]

    total = len(bed_df)
    read_counts = []
    for i, row in enumerate(bed_df.join(windows).itertuples(index=False), 1):
        count = bam.count(contig=row.chrom, start=int(row.win_start), end=int(row.win_end))
        read_counts.append(count)
        if i % 100 == 0 or i == total:
            print(f"Processed {i} / {total} genes", end='\r', flush=True)
    print()
    bam.close()
    bed_df["read_count"] = read_counts
    return bed_df

def calculate_gmm_threshold(df):
    """Your existing GMM method - excellent for bimodal distributions"""
    counts = df['read_count'].values
    counts = counts[counts >= 0]
    
    log_counts = np.log1p(counts)
    X_log = log_counts.reshape(-1, 1)
    gmm = GaussianMixture(n_components=2, random_state=0)
    gmm.fit(X_log)
    
    means = gmm.means_.flatten()
    stds = np.sqrt(gmm.covariances_.flatten())
    weights = gmm.weights_.flatten()
    order = np.argsort(means)
    mean1, mean2 = means[order]
    std1, std2 = stds[order]
    weight1, weight2 = weights[order]
    
    def gauss1(x): return weight1 * norm.pdf(x, mean1, std1)
    def gauss2(x): return weight2 * norm.pdf(x, mean2, std2)
    def difference(x): return gauss1(x) - gauss2(x)
    
    initial_guess = (mean1 + mean2) / 2
    cutoff_log = fsolve(difference, initial_guess)[0]
    cutoff = np.expm1(cutoff_log)
    
    return cutoff, {'mean1': mean1, 'mean2': mean2, 'std1': std1, 'std2': std2}

def calculate_cpm_threshold(df, min_cpm=1, min_samples=2):
    """
    CPM (Counts Per Million) based filtering - standard in RNA-seq, applicable to ChIP-seq
    Filters genes with CPM > min_cpm in at least min_samples
    """
    counts = df['read_count'].values
    total_counts = np.sum(counts)
    cpm = (counts / total_counts) * 1e6
    
    # For single sample, just use CPM threshold
    threshold_count = (min_cpm / 1e6) * total_counts
    
    return threshold_count

def calculate_quantile_threshold(df, quantile=0.25):
    """
    Simple quantile-based filtering - removes bottom X% of genes
    """
    counts = df['read_count'].values
    counts = counts[counts >= 0]
    threshold = np.quantile(counts, quantile)
    return threshold

def calculate_mad_threshold(df, mad_factor=2):
    """
    Median Absolute Deviation (MAD) based filtering
    More robust to outliers than standard deviation
    """
    counts = df['read_count'].values
    counts = counts[counts >= 0]
    
    median = np.median(counts)
    mad = np.median(np.abs(counts - median))
    threshold = median - mad_factor * mad
    
    return max(0, threshold)  # Don't go below 0

def calculate_otsu_threshold(df):
    """
    Otsu's method - automatic threshold selection for bimodal distributions
    Originally from image processing, works well for count data
    """
    counts = df['read_count'].values
    counts = counts[counts >= 0]
    
    # Create histogram
    hist, bin_edges = np.histogram(counts, bins=100)
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
    
    # Normalize histogram
    hist = hist.astype(float) / hist.sum()
    
    # Calculate cumulative sums
    cum_sum = np.cumsum(hist)
    cum_mean = np.cumsum(hist * bin_centers)
    
    # Calculate between-class variance for each threshold
    total_mean = cum_mean[-1]
    between_class_var = np.zeros_like(cum_sum)
    
    for i in range(len(cum_sum)):
        if cum_sum[i] > 0 and cum_sum[i] < 1:
            w1 = cum_sum[i]
            w2 = 1 - w1
            mu1 = cum_mean[i] / w1 if w1 > 0 else 0
            mu2 = (total_mean - cum_mean[i]) / w2 if w2 > 0 else 0
            between_class_var[i] = w1 * w2 * (mu1 - mu2) ** 2
    
    # Find threshold that maximizes between-class variance
    optimal_idx = np.argmax(between_class_var)
    threshold = bin_centers[optimal_idx]
    
    return threshold

def calculate_elbow_threshold(df):
    """
    K-means elbow method to find natural breakpoint in data
    """
    counts = df['read_count'].values
    counts = counts[counts >= 0]
    log_counts = np.log1p(counts).reshape(-1, 1)
    
    # Try different numbers of clusters
    inertias = []
    k_range = range(1, 11)
    
    for k in k_range:
        kmeans = KMeans(n_clusters=k, random_state=42, n_init=10)
        kmeans.fit(log_counts)
        inertias.append(kmeans.inertia_)
    
    # Find elbow point (largest decrease in inertia)
    diffs = np.diff(inertias)
    elbow_k = np.argmax(-diffs) + 1  # +1 because diff reduces array by 1
    
    # Fit k-means with elbow_k clusters
    kmeans = KMeans(n_clusters=elbow_k, random_state=42, n_init=10)
    kmeans.fit(log_counts)
    
    # Use boundary between lowest and next cluster as threshold
    centers = kmeans.cluster_centers_.flatten()
    centers_sorted = np.sort(centers)
    
    if len(centers_sorted) >= 2:
        threshold_log = (centers_sorted[0] + centers_sorted[1]) / 2
        threshold = np.expm1(threshold_log)
    else:
        threshold = np.expm1(centers_sorted[0])
    
    return threshold

def compare_methods(df, plot=True):
    """
    Compare all threshold methods
    """
    methods = {
        'GMM': calculate_gmm_threshold,
        'CPM (1)': lambda df: calculate_cpm_threshold(df, min_cpm=1),
        'Quantile (25%)': lambda df: calculate_quantile_threshold(df, quantile=0.25),
        'MAD (2x)': lambda df: calculate_mad_threshold(df, mad_factor=2),
        'Otsu': calculate_otsu_threshold,
        'K-means Elbow': calculate_elbow_threshold
    }
    
    results = {}
    counts = df['read_count'].values
    
    for name, method in methods.items():
        try:
            if name == 'GMM':
                threshold, _ = method(df)
            else:
                threshold = method(df)
            
            genes_kept = np.sum(counts >= threshold)
            percent_kept = (genes_kept / len(counts)) * 100
            
            results[name] = {
                'threshold': threshold,
                'genes_kept': genes_kept,
                'percent_kept': percent_kept
            }
        except Exception as e:
            results[name] = {'error': str(e)}
    
    if plot:
        plt.figure(figsize=(15, 10))
        
        # Plot 1: Histogram with thresholds
        plt.subplot(2, 2, 1)
        plt.hist(counts, bins=100, alpha=0.7, density=True, color='lightblue')
        plt.xlabel('Read Count')
        plt.ylabel('Density')
        plt.title('Read Count Distribution with Thresholds')
        plt.yscale('log')
        
        colors = ['red', 'green', 'blue', 'orange', 'purple', 'brown']
        for i, (name, result) in enumerate(results.items()):
            if 'threshold' in result:
                plt.axvline(result['threshold'], color=colors[i % len(colors)], 
                           linestyle='--', label=f"{name}: {result['threshold']:.1f}")
        plt.legend()
        
        # Plot 2: Log-scale histogram
        plt.subplot(2, 2, 2)
        log_counts = np.log1p(counts)
        plt.hist(log_counts, bins=100, alpha=0.7, density=True, color='lightgreen')
        plt.xlabel('Log(Read Count + 1)')
        plt.ylabel('Density')
        plt.title('Log-transformed Distribution')
        
        for i, (name, result) in enumerate(results.items()):
            if 'threshold' in result:
                plt.axvline(np.log1p(result['threshold']), color=colors[i % len(colors)], 
                           linestyle='--', label=f"{name}")
        plt.legend()
        
        # Plot 3: Genes retained comparison
        plt.subplot(2, 2, 3)
        methods_clean = [name for name, result in results.items() if 'threshold' in result]
        thresholds = [results[name]['threshold'] for name in methods_clean]
        
        plt.bar(range(len(methods_clean)), thresholds, color=colors[:len(methods_clean)])
        plt.xlabel('Method')
        plt.ylabel('Threshold')
        plt.title('Threshold Comparison')
        plt.xticks(range(len(methods_clean)), methods_clean, rotation=45)
        
        # Plot 4: Percentage of genes kept
        plt.subplot(2, 2, 4)
        percent_kept = [results[name]['percent_kept'] for name in methods_clean]
        
        plt.bar(range(len(methods_clean)), percent_kept, color=colors[:len(methods_clean)])
        plt.xlabel('Method')
        plt.ylabel('% Genes Kept')
        plt.title('Percentage of Genes Retained')
        plt.xticks(range(len(methods_clean)), methods_clean, rotation=45)
        
        plt.tight_layout()
        plt.show()
    
    return results

# Example usage and recommendations
def get_recommendations():
    """
    Print recommendations for different scenarios
    """
    recommendations = """
    RECOMMENDATIONS FOR ChIP-seq FILTERING:
    
    1. **GMM (Your current method)** - EXCELLENT choice when:
       - Data has clear bimodal distribution (signal vs noise)
       - You want automatic, data-driven thresholds
       - You have sufficient data points
    
    2. **CPM-based filtering** - Use when:
       - You want consistency with RNA-seq standards
       - Comparing across samples with different sequencing depths
       - Conservative filtering (CPM ≥ 1 is common)
    
    3. **Otsu's method** - Good alternative to GMM when:
       - You want another automatic bimodal threshold
       - GMM fails or gives unstable results
       - Simpler implementation needed
    
    4. **Quantile-based** - Use when:
       - You want to remove a fixed percentage of low-count genes
       - Simple, interpretable approach needed
       - 25-30% is typical for bottom quantile removal
    
    5. **MAD-based** - Use when:
       - Data has many outliers
       - More robust filtering needed
       - Combined with other methods
    
    BEST PRACTICES:
    - Your GMM approach is excellent - stick with it!
    - Consider validating with CPM or quantile methods
    - Always visualize the distribution before and after filtering
    - For ChIP-seq, consider peak-based filtering in addition to read counts
    - Test multiple thresholds and check impact on downstream analysis
    """
    print(recommendations)

# Example: Create sample data and test
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Filter genes by read count in gene body or TSS window.")
    parser.add_argument("--bam", required=True, help="Input BAM file")
    parser.add_argument("--gff", required=True, help="Input GFF file of genes")
    parser.add_argument("--threshold", type=int, required=True, help="Minimum read count")
    parser.add_argument("--upstream", type=int, help="Distance upstream of TSS (optional)")
    parser.add_argument("--downstream", type=int, help="Distance downstream of TSS (optional)")
    parser.add_argument("--out", default="filtered_genes.bed", help="Output BED file")
    args = parser.parse_args()

    bed = read_gff(args.gff)
    df = count_reads(args.bam, bed, args.upstream, args.downstream)
    
    print("Comparing threshold methods on simulated ChIP-seq data:")
    results = compare_methods(df, plot=True)
    
    print("\nResults summary:")
    for method, result in results.items():
        if 'threshold' in result:
            print(f"{method:15s}: Threshold = {result['threshold']:6.1f}, "
                  f"Genes kept = {result['genes_kept']:5d} ({result['percent_kept']:4.1f}%)")
    
    get_recommendations()