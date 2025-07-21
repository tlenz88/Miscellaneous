import cooler
import cooltools
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from scipy import stats
import pandas as pd

def parse_region(region_str):
    """
    Parse region string in format 'chr1:15000000-20000000'
    Returns tuple of (chrom, start, end)
    """
    chrom = region_str.split(':')[0]
    start = int(region_str.split(':')[1].split('-')[0])
    end = int(region_str.split(':')[1].split('-')[1])
    return chrom, start, end

def load_and_process_cool(cool_file, region_str):
    """
    Load and extract data from .cool file for a specific region
    """
    # Parse region
    chrom, start, end = parse_region(region_str)
    
    # Load the cooler file
    c = cooler.Cooler(cool_file)
    
    # Get matrix for specified region
    matrix = c.matrix(balance=True).fetch(f"{chrom}:{start}-{end}")
    
    # Replace NaN values with 0
    matrix = np.nan_to_num(matrix, nan=0.0)
    
    return matrix

def calculate_oe_matrix(matrix):
    """
    Calculate observed/expected matrix
    """
    n = len(matrix)
    expected = np.zeros(n)
    for i in range(n):
        diag = np.diag(matrix, k=i)
        expected[i] = np.mean(diag[diag != 0])
    
    expected_matrix = np.zeros_like(matrix)
    for i in range(n):
        for j in range(n):
            expected_matrix[i,j] = expected[abs(i-j)]
    
    oe_matrix = np.divide(matrix, expected_matrix, 
                         out=np.zeros_like(matrix), 
                         where=expected_matrix!=0)
    
    return oe_matrix

def generate_correlation_matrix(cool_file, region_str):
    """
    Generate correlation matrix using cooltools methods
    """
    # Parse region
    chrom, start, end = parse_region(region_str)
    
    clr = cooler.Cooler(cool_file)
    
    # Create view_df for the specific region
    view_df = pd.DataFrame({
        'chrom': [chrom],
        'start': [start],
        'end': [end],
        'name': [f"{chrom}_{start}_{end}"]
    })

    # Get expected values
    expected_df = cooltools.expected_cis(
        clr,
        view_df=view_df,
        clr_weight_name="weight",
        nproc=8
    )
    
    # Get observed matrix for the region
    observed = clr.matrix(balance=True).fetch(f"{chrom}:{start}-{end}")
    observed = np.nan_to_num(observed, nan=0.0)
    
    n = observed.shape[0]
    expected_matrix = np.zeros_like(observed)
    
    # Fill in expected matrix using smoothed values
    for idx, row in expected_df.iterrows():
        diag_idx = int(row['dist']/clr.binsize)
        if diag_idx < n:
            exp_val = row['balanced.avg.smoothed']
            if pd.isna(exp_val):
                exp_val = row['balanced.avg.smoothed.agg']
            if pd.isna(exp_val):
                continue
                
            for i in range(n-diag_idx):
                expected_matrix[i, i+diag_idx] = exp_val
                if diag_idx > 0:
                    expected_matrix[i+diag_idx, i] = exp_val
    
    # Handle zeros in expected matrix
    expected_matrix[expected_matrix == 0] = np.min(expected_matrix[expected_matrix > 0])
    
    # Calculate O/E matrix
    oe_matrix = observed / expected_matrix
    
    # Handle extreme values
    oe_matrix = np.clip(oe_matrix, 0, np.percentile(oe_matrix[~np.isnan(oe_matrix) & ~np.isinf(oe_matrix)], 95))
    
    # Replace remaining inf and nan values with 0
    oe_matrix = np.nan_to_num(oe_matrix, nan=0.0, posinf=0.0, neginf=0.0)
    
    # Calculate correlation
    corr_matrix = np.corrcoef(oe_matrix)
    
    return corr_matrix


def plot_correlation_heatmap(matrix, title, ax):
    """
    Plot correlation heatmap using seaborn with symmetric log scale
    """
    if np.all(np.isnan(matrix)):
        ax.text(0.5, 0.5, 'Invalid matrix', ha='center', va='center')
        return
    
    # Create a symmetric log norm
    from matplotlib.colors import SymLogNorm
    
    # Define the normalization
    linthresh = 0.05
    norm = SymLogNorm(linthresh=linthresh, vmin=-0.5, vmax=0.5)
    
    # Create heatmap with symlog scale
    sns.heatmap(matrix, 
                cmap='RdBu_r',
                norm=norm,
                center=0,
                ax=ax,
                rasterized=True,
                cbar_kws={'label': 'Correlation (symlog scale)'})
    
    ax.set_title(title)
    
def main(cool_file, region_str, output_pdf):
    # Load data for region
    matrix = load_and_process_cool(cool_file, region_str)
    
    # Calculate correlations
    cooltools_corr = generate_correlation_matrix(cool_file, region_str)
    print(np.min(cooltools_corr))
    print(np.max(cooltools_corr))
    oe_matrix = calculate_oe_matrix(matrix)
    manual_corr = np.corrcoef(oe_matrix)
    print(np.min(manual_corr))
    print(np.max(manual_corr))
    direct_corr = np.corrcoef(matrix)
    print(np.min(direct_corr))
    print(np.max(direct_corr))

    # Create figure with more space for colorbars
    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(18, 5), dpi=300)
    
    # Plot heatmaps
    plot_correlation_heatmap(cooltools_corr, 'Cooltools Correlation', ax1)
    plot_correlation_heatmap(manual_corr, 'Manual O/E Correlation', ax2)
    plot_correlation_heatmap(direct_corr, 'Direct Correlation', ax3)
    
    plt.tight_layout()
    
    # Save figure
    plt.savefig(output_pdf, 
                format='pdf', 
                bbox_inches='tight', 
                dpi=300)
    plt.close()

if __name__ == "__main__":
    import sys
    if len(sys.argv) != 4:
        print("Usage: python script.py cool_file region output_pdf")
        print("Region format example: chr1:15000000-20000000")
        sys.exit(1)
    
    cool_file = sys.argv[1]
    region_str = sys.argv[2]
    output_pdf = sys.argv[3]
    main(cool_file, region_str, output_pdf)