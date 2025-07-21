import bioframe
import cooler
import cooler.balance
import cooltools
import sys
import warnings
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats


# Filter out RuntimeWarnings from numpy
warnings.filterwarnings('ignore', category=RuntimeWarning)

def assign_colors(merged_df):
    """
    Assign colors based on distance from diagonal and compartment type.
    """
    # Ensure no NaN values
    df = merged_df.dropna(subset=['E1z_df1', 'E1z_df2']).copy()
    
    # Calculate differences and standard deviation
    df['z_diff'] = df['E1z_df1'] - df['E1z_df2']
    diagonal_dist = 2 * np.nanstd(df['z_diff'])
    
    # Initialize colors array
    colors = np.full(len(df), 'gray', dtype=object)
    
    # A -> B (light blue): starts positive, significant decrease
    a_to_b = (df['E1z_df1'] > 0) & (df['E1z_df1'] > df['E1z_df2'] + diagonal_dist) & (df['E1z_df2'] < 0)
    colors[a_to_b] = '#56B4E9'
    
    # A -> weaker A (orange): starts positive, moderate decrease
    a_to_weaker = (df['E1z_df1'] > 0) & (df['E1z_df1'] > df['E1z_df2'] + diagonal_dist) & (df['E1z_df2'] > 0)
    colors[a_to_weaker] = '#E69F00'
    
    # A -> stronger A (green): starts positive, significant increase
    a_to_stronger = (df['E1z_df1'] > 0) & (df['E1z_df2'] > df['E1z_df1'] + diagonal_dist) & (df['E1z_df2'] > 0)
    colors[a_to_stronger] = '#009E73'
    
    # B -> A (blue): starts negative, significant increase
    b_to_a = (df['E1z_df1'] < 0) & (df['E1z_df2'] > df['E1z_df1'] + diagonal_dist) & (df['E1z_df2'] > 0)
    colors[b_to_a] = '#0072B2'
    
    # B -> weaker B (red): starts negative, significant decrease
    b_to_weaker = (df['E1z_df1'] < 0) & (df['E1z_df1'] < df['E1z_df2'] - diagonal_dist) & (df['E1z_df2'] < 0)
    colors[b_to_weaker] = '#D55E00'
    
    # B -> stronger B (pink): starts negative, moderate increase
    b_to_stronger = (df['E1z_df1'] < 0) & (df['E1z_df2'] < df['E1z_df1'] - diagonal_dist) & (df['E1z_df2'] < 0)
    colors[b_to_stronger] = '#CC79A7'
    
    return colors


def calculate_compartments(clr: cooler.Cooler, view_df) -> pd.DataFrame:
    """
    Perform compartment analysis.
    """
    print("Calculating compartments...")
    bins = clr.bins()[:]
    hg38_genome = bioframe.load_fasta('/mnt/f/organism_genome/Hsapien_autosomes/GRCh38.fasta')
    gc_track = bioframe.frac_gc(bins[['chrom', 'start', 'end']], hg38_genome)
    eigvals, eigvecs = cooltools.eigs_cis(
        clr,
        view_df=view_df,
        n_eigs=2,
        phasing_track=gc_track,
        ignore_diags=2,
        sort_metric=None
    )
    eigs_df = pd.DataFrame({
        'chrom': np.repeat(clr.chromnames, clr.bins()[:].groupby('chrom').size()),
        'start': bins['start'],
        'end': bins['end'],
        'E1': eigvecs['E1'],
        'E2': eigvecs['E2']
    })
    return eigs_df

def balance_matrix(clr: cooler.Cooler) -> cooler.Cooler:
    """
    Perform matrix balancing if not already balanced.
    """
    if not clr.info['metadata'].get('balance-applied', False):
        cooler.balance_cooler(clr)
    return clr

def plot_compartment_comparison(file1_path, file2_path, output_path, zscore_cutoff=4):
    """
    Create a scatter plot comparing compartment scores between conditions.
    """
    clr1 = cooler.Cooler(file1_path)
    clr_balanced1 = balance_matrix(clr1)
    view_df1 = pd.DataFrame({'chrom': clr_balanced1.chromnames,
                            'start': 0,
                            'end': clr_balanced1.chromsizes.values,
                            'name': clr_balanced1.chromnames}
                           )
    eigs_df1 = calculate_compartments(clr_balanced1, view_df1)
    eigs_df1 = eigs_df1.dropna(subset=['E1'])

    clr2 = cooler.Cooler(file2_path)
    clr_balanced2 = balance_matrix(clr2)
    view_df2 = pd.DataFrame({'chrom': clr_balanced2.chromnames,
                            'start': 0,
                            'end': clr_balanced2.chromsizes.values,
                            'name': clr_balanced2.chromnames}
                           )
    eigs_df2 = calculate_compartments(clr_balanced2, view_df2)
    eigs_df2 = eigs_df2.dropna(subset=['E1'])
    
    merged_df = pd.merge(eigs_df1[['chrom', 'start', 'end', 'E1']],
                         eigs_df2[['chrom', 'start', 'end', 'E1']],
                         on=['chrom', 'start', 'end'],
                         suffixes=('_df1', '_df2'))
    
    merged_df['E1z_df1'] = stats.zscore(merged_df['E1_df1'])
    merged_df['E1z_df2'] = stats.zscore(merged_df['E1_df2'])
    
    # Assign colors
    colors = assign_colors(merged_df)
    
    # Calculate correlation
    correlation = np.corrcoef(merged_df['E1z_df1'], merged_df['E1z_df2'])[0, 1]
    
    # Calculate plot limits
    max_abs_value = max(
        np.nanmax(np.abs(merged_df['E1z_df1'])),
        np.nanmax(np.abs(merged_df['E1z_df2']))
    )
    
    # Create plot
    plt.figure(figsize=(10, 10))
    
    # Add vertical and horizontal lines at 0
    plt.axhline(y=0, color='black', linestyle='--')
    plt.axvline(x=0, color='black', linestyle='--')
    
    # Create scatter plot with colors
    plt.scatter(merged_df['E1z_df1'], merged_df['E1z_df2'], c=colors, s=5)
    
    # Set equal limits centered at 0
    plt.xlim(-max_abs_value, max_abs_value)
    plt.ylim(-max_abs_value, max_abs_value)
    
    # Add labels and title
    plt.xlabel('Control compartment score (E1)')
    plt.ylabel('Treatment compartment score (E1)')
    plt.title(f'Compartment change\nr = {correlation:.3f}')
    
    # Calculate category counts and percentages
    unique_colors = {
        '#56B4E9': 'A -> B',
        '#E69F00': 'A -> weaker A',
        '#009E73': 'A -> stronger A',
        '#0072B2': 'B -> A',
        '#D55E00': 'B -> weaker B',
        '#CC79A7': 'B -> stronger B'
    }
    total_points = len(colors)
    category_counts = {color: (colors == color).sum() for color in unique_colors.keys()}
    category_percentages = {
        color: count / total_points * 100 for color, count in category_counts.items()
    }
    
    # Create legend with percentages
    legend_elements = [
        plt.Line2D(
            [0], [0], marker='o', color='w', markerfacecolor=color,
            label=f"{label} ({category_percentages[color]:.1f}%)", markersize=8
        )
        for color, label in unique_colors.items()
    ]
    plt.legend(handles=legend_elements, loc='center left', bbox_to_anchor=(1, 0.5))
    
    # Ensure square aspect ratio
    plt.axis('square')
    
    # Save plot with extra space for legend
    plt.savefig(output_path, bbox_inches='tight')
    plt.close()

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python script.py <control.cool> <treatment.cool> <output.pdf>")
        sys.exit(1)
    try:
        plot_compartment_comparison(sys.argv[1], sys.argv[2], sys.argv[3])
    except Exception as e:
        print(f"Error: {str(e)}")
        sys.exit(1)
