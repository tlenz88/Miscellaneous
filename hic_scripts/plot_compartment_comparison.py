import sys
import numpy as np
import matplotlib.pyplot as plt
import pyBigWig
from scipy import stats
from scipy.stats import wilcoxon

def get_interval_values(bw, chrom, window_size=100000):
    """
    Get values from bigWig intervals directly, binning where needed.
    """
    intervals = bw.intervals(chrom)
    if not intervals:
        return np.array([])
    
    chrom_length = bw.chroms()[chrom]
    n_bins = int(np.ceil(chrom_length / window_size))
    binned_values = np.zeros(n_bins)
    
    for start, end, value in intervals:
        bin_idx = start // window_size
        if bin_idx < n_bins:
            binned_values[bin_idx] = value
    
    return binned_values

def assign_colors(values1, values2):
    """
    Assign colors based on distance from diagonal and compartment type.
    
    Parameters:
    -----------
    values1, values2 : numpy.ndarray
        Arrays of compartment scores
    """
    # Calculate distance from diagonal
    diagonal_dist = 1.5 * np.std(values2 - values1)

    # Initialize colors array
    colors = ['gray'] * len(values1)
    
    # A -> B (light blue): starts positive, significant decrease
    a_to_b = (values1 > 0) & (values1 > values2 + diagonal_dist) & (values2 < 0)
    for idx in np.where(a_to_b)[0]:
        colors[idx] = '#56B4E9'
    
    # A -> weaker A (orange): starts positive, moderate decrease
    a_to_weaker = (values1 > 0) & (values1 > values2 + diagonal_dist) & (values2 > 0)
    for idx in np.where(a_to_weaker)[0]:
        colors[idx] = '#E69F00'
    
    # A -> stronger A (green): starts positive, significant increase
    a_to_stronger = (values1 > 0) & (values2 > values1 + diagonal_dist) & (values2 > 0)
    for idx in np.where(a_to_stronger)[0]:
        colors[idx] = '#009E73'
    
    # B -> A (blue): starts negative, significant increase
    b_to_a = (values1 < 0) & (values2 > values1 + diagonal_dist) & (values2 > 0)
    for idx in np.where(b_to_a)[0]:
        colors[idx] = '#0072B2'

    # B -> weaker B (red): starts negative, significant decrease
    b_to_weaker = (values1 < 0) & (values1 < values2 - diagonal_dist) & (values2 < 0)
    for idx in np.where(b_to_weaker)[0]:
        colors[idx] = '#D55E00'

    # B -> stronger B (pink): starts negative, moderate increase
    b_to_stronger = (values1 < 0) & (values2 < values1 - diagonal_dist) & (values2 < 0)
    for idx in np.where(b_to_stronger)[0]:
        colors[idx] = '#CC79A7'

    return colors

def plot_compartment_comparison(bigwig1_path, bigwig2_path, output_path, 
                                window_size=100000, zscore_cutoff=4):
    """
    Create a scatter plot comparing compartment scores between conditions.
    """
    bw1 = pyBigWig.open(bigwig1_path)
    bw2 = pyBigWig.open(bigwig2_path)
    
    values1 = []
    values2 = []
    
    # Process each chromosome
    for chrom in bw1.chroms().keys():
        print(f"Processing {chrom}...")
        
        binned1 = get_interval_values(bw1, chrom, window_size)
        binned2 = get_interval_values(bw2, chrom, window_size)
        
        min_length = min(len(binned1), len(binned2))
        binned1 = binned1[:min_length]
        binned2 = binned2[:min_length]
        
        non_zero_mask = (binned1 != 0) | (binned2 != 0)
        values1.extend(binned1[non_zero_mask])
        values2.extend(binned2[non_zero_mask])
    
    values1 = np.array(values1)
    values2 = np.array(values2)
    
    # Remove outliers
    z1 = np.abs(stats.zscore(values1))
    z2 = np.abs(stats.zscore(values2))
    valid_points = (z1 < zscore_cutoff) & (z2 < zscore_cutoff)
    values1 = values1[valid_points]
    values2 = values2[valid_points]
    
    # Assign colors based on distance from diagonal
    colors = assign_colors(values1, values2)
    
    # Calculate correlation
    correlation = np.corrcoef(values1, values2)[0, 1]
    
    # Calculate plot limits
    max_abs_value = max(
        max(values1.max(), abs(values1.min())),
        max(values2.max(), abs(values2.min())))
    
    # Create plot
    plt.figure(figsize=(10, 10))
    
    # Add vertical and horizontal lines at 0
    plt.axhline(y=0, color='black', linestyle='--')
    plt.axvline(x=0, color='black', linestyle='--')
    
    # Create scatter plot with colors
    plt.scatter(values1, values2, c=colors, s=5)
    
    # Set equal limits centered at 0
    plt.xlim(-max_abs_value, max_abs_value)
    plt.ylim(-max_abs_value, max_abs_value)
    
    # Add labels and title
    plt.xlabel('HFF compartment score (E1)')
    plt.ylabel('ME49RFP_24hpi compartment score (E1)')
    window_size_kb = window_size / 1000
    plt.title(f'Compartment change\nr = {correlation:.3f}')
    
    # Add legend
    legend_elements = [
        plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='#56B4E9', label='A -> B', markersize=8),
        plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='#E69F00', label='A -> weaker A', markersize=8),
        plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='#009E73', label='A -> stronger A', markersize=8),
        plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='#0072B2', label='B -> A', markersize=8),
        plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='#D55E00', label='B -> weaker B', markersize=8),
        plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='#CC79A7', label='B -> stronger B', markersize=8)
    ]
    plt.legend(handles=legend_elements, loc='center left', bbox_to_anchor=(1, 0.5))

    # Ensure square aspect ratio
    plt.axis('square')
    
    # Save plot with extra space for legend
    plt.savefig(output_path, bbox_inches='tight')
    plt.close()
    
    # Close bigWig files
    bw1.close()
    bw2.close()
    
    # Print statistics
    print("\nCompartment change statistics:")
    categories = {
        'A->B': np.sum(colors == 'red'),
        'A->weaker A': np.sum(colors == 'salmon'),
        'B->stronger B': np.sum(colors == 'orange'),
        'B->A': np.sum(colors == 'darkblue'),
        'B->weaker B': np.sum(colors == 'lightblue'),
        'A->stronger A': np.sum(colors == 'dodgerblue'),
        'No significant change': np.sum(colors == 'gray')
    }
    
    total_points = len(colors)
    for category, count in categories.items():
        percentage = (count / total_points) * 100
        print(f"{category}: {count} regions ({percentage:.1f}%)")
    
    return values1, values2

if __name__ == "__main__":
    # First two args should be bigwig files output by hicPCA
    plot_compartment_comparison(sys.argv[1], sys.argv[2], f'ME49RFP_24hpi_vs_HFF_{sys.argv[3]}_compartment_comparison.pdf', window_size=int(sys.argv[3]))