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
from matplotlib import ticker
from scipy.stats import pearsonr
import cooltools.lib.plotting
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.backends.backend_pdf import PdfPages

def load_data(expected_file, eigs_file, cool_file):
    """Load and process input files"""
    expected_df = pd.read_csv(expected_file, sep='\t')
    eigs_df = pd.read_csv(eigs_file, sep='\t')
    clr = cooler.Cooler(cool_file)
    
    hg38_chromsizes = bioframe.fetch_chromsizes('hg38')
    hg38_cens = bioframe.fetch_centromeres('hg38')
    hg38_arms = bioframe.make_chromarms(hg38_chromsizes,  hg38_cens)
    hg38_arms = hg38_arms[hg38_arms.chrom.isin(clr.chromnames)].reset_index(drop=True)

    # Convert base pairs to Mb for better plotting
    expected_df['dist_mb'] = expected_df['dist_bp'] / 1e6
    
    return expected_df, eigs_df, clr, hg38_arms

def get_compartment_track(eigs_df):
    """Create compartment track for visualization"""
    track = eigs_df[['chrom', 'start', 'end', 'E1']].copy()
    track.loc[:, 'compartment'] = np.where(track['E1'] > 0, 1, -1)
    return track

def plot_heatmaps(clr, chrom, eigs_df, resolution):
    chrom_length = clr.chromsizes[chrom]
    track = get_compartment_track(eigs_df)
    
    mat = clr.matrix(balance=True).fetch(chrom)
    n = mat.shape[0]

    # Identify rows and columns where all entries are zero
    non_zero_rows = np.any(mat != 0, axis=1)
    non_zero_cols = np.any(mat != 0, axis=0)

    # Filter out rows and columns where all entries are zero (removes rows/columns with all zeros)
    mat_filtered = mat[non_zero_rows, :][:, non_zero_cols]

    # Calculate the expected matrix for the filtered matrix
    expected = np.zeros_like(mat_filtered)
    for i in range(mat_filtered.shape[0]):
        diag_values = np.diagonal(mat_filtered, i)
        valid_values = diag_values[diag_values > 0]
        if len(valid_values) > 0:
            diag_mean = np.mean(valid_values)
            indices = np.indices(mat_filtered.shape)
            expected[indices[0], indices[1]] = np.where(
                np.abs(indices[0] - indices[1]) == i,
                diag_mean,
                expected[indices[0], indices[1]]
            )

    # Calculate observed/expected matrix for the filtered matrix
    obs_exp = np.zeros_like(mat_filtered)
    mask = (expected > 0) & (mat_filtered > 0)
    obs_exp[mask] = mat_filtered[mask] / expected[mask]
    obs_exp = np.nan_to_num(obs_exp, nan=0, posinf=0, neginf=0)
    
    obs_exp_log = np.log10(obs_exp, where=obs_exp>0)

    fig = plt.figure(figsize=(10, 12))
    gs = gridspec.GridSpec(2, 1, height_ratios=[4, 1], hspace=0.3)
    
    ax_heatmap = fig.add_subplot(gs[0])
    start, end = 0, mat_filtered.shape[0]
    pos_window = np.linspace(0, chrom_length, mat_filtered.shape[0])[start:end]

    min_val = np.nanmin(obs_exp_log)
    max_val = np.nanmax(obs_exp_log)
    vminmax = max(-min_val, max_val)
    norm = colors.SymLogNorm(linthresh=.5, linscale=1, vmin=-vminmax, vmax=vminmax, base=10)

    im = ax_heatmap.imshow(obs_exp_log[start:end, start:end], cmap='RdBu_r', 
                           norm=norm,
                           extent=[pos_window[0]/1e6, pos_window[-1]/1e6, 
                                   pos_window[-1]/1e6, pos_window[0]/1e6])
    cbar = fig.colorbar(im, ax=ax_heatmap, label='Observed/Expected')
    ax_heatmap.set_title(f"Obs/Exp Heatmap for {chrom}")
    ax_heatmap.set_ylabel("Position (Mb)")

    # Compartment track plot
    track_minmax = max(-np.nanmin(track['E1']), np.nanmax(track['E1']))
    ax_track = fig.add_subplot(gs[1], sharex=ax_heatmap)
    ax_track.plot(track['start'] / 1e6, track['E1'], color='black', lw=0)
    ax_track.fill_between(track['start'] / 1e6, 0, track['E1'],
                          where=track['E1'] > 0, facecolor='red', alpha=0.5)
    ax_track.fill_between(track['start'] / 1e6, 0, track['E1'],
                          where=track['E1'] < 0, facecolor='blue', alpha=0.5)
    ax_track.set_ylim(-track_minmax, track_minmax)
    #ax_track.set_yticks([-1, 1])
    #ax_track.set_yticklabels(['-1', '1'])
    ax_track.set_xlabel("Position (Mb)")
    ax_track.set_ylabel("Compartment")

    return fig

def main():
    parser = argparse.ArgumentParser(description='Analyze Hi-C compartments and scaling')
    parser.add_argument('--expected', required=True, help='Expected contacts TSV file')
    parser.add_argument('--eigs', required=True, help='Eigenvectors TSV file')
    parser.add_argument('--cool', required=True, help='Balanced .cool file')
    parser.add_argument('--output', required=True, help='Output directory')
    parser.add_argument('--resolution', required=False, help='Binning resolution')
    args = parser.parse_args()

    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)
    resolution = int(args.resolution)

    expected_df, eigs_df, clr, hg38_arms = load_data(args.expected, args.eigs, args.cool)
    outfile = output_dir / f'obs_exp_{resolution}_heatmap.pdf'
    with PdfPages(outfile) as pdf:
        for chrom in eigs_df['chrom'].unique():
            chrom_eigs = eigs_df[eigs_df['chrom'] == chrom]
       	    #fig = plot_heatmaps(clr, chrom, chrom_eigs, resolution)
            fig = plot_heatmaps(clr, chrom, expected_df, resolution)
       	    pdf.savefig(fig)
       	    plt.close(fig)

if __name__ == "__main__":
    main()