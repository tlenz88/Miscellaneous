import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import argparse
from pathlib import Path
import cooler
import cooltools
import bioframe
import itertools
import matplotlib.colors as colors
from matplotlib.colors import Normalize, LogNorm
import matplotlib.gridspec as gridspec
from matplotlib.gridspec import GridSpec
import cooltools.lib.plotting
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.backends.backend_pdf import PdfPages

def find_region_eigs(sample_eigs, region):
    region_eigs = pd.read_csv(sample_eigs, sep='\t')
    region_eigs.loc[np.isnan(region_eigs['weight']), 'weight'] = 0
    region_eigs.loc[np.isnan(region_eigs['E1']), 'E1'] = 0

    try:
        chrom, st_end = region.split(':')
        start, end = st_end.split('-')
        region_eigs_filtered = region_eigs[(region_eigs['chrom'] == chrom) & 
        (region_eigs['start'] < int(end)) & 
        (region_eigs['end'] > int(start))]
    except:
        chrom = region
        region_eigs_filtered = region_eigs[(region_eigs['chrom'] == chrom)]

    return region_eigs_filtered

def plot_compartments(eigs_df):
    """Plot compartment eigenvectors"""
    fig, axes = plt.subplots(len(eigs_df), 1, figsize=(15, 10), sharex=True)
    
    if len(eigs_df) == 1:
        axes = [axes]

    max_E1 = max(df['E1'].max() for df in eigs_df.values())

    for ax, (key, df) in zip(axes, eigs_df.items()):
        pos = df['start'] / 1e6
        ax.plot(df['start'] / 1e6, df['E1'], linewidth=0.5)
        ax.set_ylabel(key)
        ax.grid(True, alpha=0.3)
        ax.set_ylim(-max_E1, max_E1)
        ax.fill_between(df['start'] / 1e6, 0, df['E1'], where=(df['E1'] > 0), color='red', alpha=0.5)
        ax.fill_between(df['start'] / 1e6, 0, df['E1'], where=(df['E1'] < 0), color='blue', alpha=0.5)
        ax.axhline(0, color='black', linewidth=0.5)

    axes[-1].set_xlabel('Genomic Position (Mb)')
    plt.suptitle(f'Compartment Analysis')
    
    return fig

def main():
    parser = argparse.ArgumentParser(description='Analyze Hi-C compartments and scaling')
    parser.add_argument('--eigs', required=True, nargs='+', help='Eigenvectors TSV file')
    parser.add_argument('--output', required=True, help='Output file')
    parser.add_argument('--region', required=True, help='Genome region')
    parser.add_argument('--samples', nargs='+', required=False)
    args = parser.parse_args()

    eigs_df = {}
    for idx, sample in enumerate(args.samples):
        eigs_df[sample] = find_region_eigs(args.eigs[idx], args.region)

    with PdfPages(args.output) as pdf:
        track_fig = plot_compartments(eigs_df)
        pdf.savefig(track_fig)
        plt.close(track_fig)

if __name__ == "__main__":
    main()