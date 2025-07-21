import bioframe
import cooler
import cooltools
import itertools
import sys
import cooltools.lib.plotting
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from cooler import balance
from cooltools import insulation
from matplotlib.colors import LogNorm
from matplotlib.ticker import EngFormatter
from mpl_toolkits.axes_grid1 import make_axes_locatable
from pathlib import Path
from scipy.stats import describe
from matplotlib.backends.backend_pdf import PdfPages

def balance_matrix(clr: cooler.Cooler) -> cooler.Cooler:
    if not clr.info['metadata'].get('balance-applied', False):
        cooler.balance_cooler(clr)
    return clr

# Initialize cooler file and metadata
clr_balanced = cooler.Cooler(sys.argv[1])
#clr_balanced = balance_matrix(clr)
resolution = clr_balanced.binsize
sample = Path(sys.argv[1]).stem

windows = [3 * resolution, 5 * resolution, 10 * resolution, 25 * resolution]

# Calculate insulation scores
insulation_table = insulation(clr_balanced, windows)
"""
# Initialize list to store chromosome-wise statistics
stats_list = []

# Initialize lists to collect all insulation scores and boundary strengths for genome-wide summary
all_insulation_scores = []
all_boundary_strengths = []
all_boundary_distances = []

# Iterate over all chromosomes in the cooler file
chromsizes = clr_balanced.chromsizes
for chrom in chromsizes.index:
    chrom_size = chromsizes[chrom]
    region = (chrom, 0, chrom_size)

    # Extract insulation data for the chromosome
    insul_region = bioframe.select(insulation_table, region)

    # Exclude NaN values for insulation scores and boundary strengths
    valid_scores = insul_region[f'log2_insulation_score_{windows[0]}'].dropna()
    valid_boundaries = insul_region[~np.isnan(insul_region[f'boundary_strength_{windows[0]}'])]

    # Calculate statistics for insulation scores
    score_stats = describe(valid_scores)
    boundary_strengths = valid_boundaries[f'boundary_strength_{windows[0]}']
    boundary_positions = valid_boundaries[['start', 'end']].mean(axis=1)

    # Calculate mean and median distance between boundaries
    mean_distance = np.mean(np.diff(boundary_positions)) if len(boundary_positions) > 1 else np.nan
    median_distance = np.median(np.diff(boundary_positions)) if len(boundary_positions) > 1 else np.nan

    # Append the chromosome-wide statistics to the list
    stats_list.append({
        'Chromosome': chrom,
        'Mean Insulation Score': score_stats.mean,
        'Median Insulation Score': np.median(valid_scores),
        '25th Percentile Insulation Score': np.percentile(valid_scores, 25),
        '75th Percentile Insulation Score': np.percentile(valid_scores, 75),
        'Min Insulation Score': score_stats.minmax[0],
        'Max Insulation Score': score_stats.minmax[1],
        'Mean Boundary Strength': np.mean(boundary_strengths),
        'Median Boundary Strength': np.median(boundary_strengths),
        'Mean Distance Between Boundaries': mean_distance,
        'Median Distance Between Boundaries': median_distance
    })

    # Collect data for genome-wide summary
    all_insulation_scores.extend(valid_scores)
    all_boundary_strengths.extend(boundary_strengths)
    all_boundary_distances.extend(np.diff(boundary_positions))

# Convert the list of stats into a DataFrame
stats_df = pd.DataFrame(stats_list)

# Calculate genome-wide summary statistics
genomewide_stats = {
    'Chromosome': 'Genome-wide',
    'Mean Insulation Score': np.mean(all_insulation_scores),
    'Median Insulation Score': np.median(all_insulation_scores),
    '25th Percentile Insulation Score': np.percentile(all_insulation_scores, 25),
    '75th Percentile Insulation Score': np.percentile(all_insulation_scores, 75),
    'Min Insulation Score': np.min(all_insulation_scores),
    'Max Insulation Score': np.max(all_insulation_scores),
    'Mean Boundary Strength': np.mean(all_boundary_strengths),
    'Median Boundary Strength': np.median(all_boundary_strengths),
    'Mean Distance Between Boundaries': np.mean(all_boundary_distances) if all_boundary_distances else np.nan,
    'Median Distance Between Boundaries': np.median(all_boundary_distances) if all_boundary_distances else np.nan
}

# Convert the genome-wide summary into a DataFrame and append it to the main stats DataFrame
genomewide_df = pd.DataFrame([genomewide_stats])
final_stats_df = pd.concat([stats_df, genomewide_df], ignore_index=True)

# Save statistics to a CSV file
output_csv = f'{Path(sys.argv[1]).parent.resolve()}/{sample}_insulation_statistics.csv'
final_stats_df.to_csv(output_csv, index=False)

print(f"Genome-wide insulation statistics saved to: {output_csv}")
"""
def pcolormesh_45deg(ax, matrix_c, start=0, resolution=1, *args, **kwargs):
    start_pos_vector = [start+resolution*i for i in range(len(matrix_c)+1)]

    n = matrix_c.shape[0]
    t = np.array([[1, 0.5], [-1, 0.5]])
    matrix_a = np.dot(np.array([(i[1], i[0])
                                for i in itertools.product(start_pos_vector[::-1],
                                                           start_pos_vector)]), t)
    x = matrix_a[:, 1].reshape(n + 1, n + 1)
    y = matrix_a[:, 0].reshape(n + 1, n + 1)
    im = ax.pcolormesh(x, y, np.flipud(matrix_c), *args, **kwargs)
    im.set_rasterized(True)
    return im


def format_ticks(ax, x=True, y=True, rotate=True):
    if y:
        ax.yaxis.set_major_formatter(bp_formatter)
    if x:
        ax.xaxis.set_major_formatter(bp_formatter)
        ax.xaxis.tick_bottom()
    if rotate:
        ax.tick_params(axis='x',rotation=45)


bp_formatter = EngFormatter('b')
plt.rcParams['font.size'] = 12
"""
regions = {'IL6': ['chr7', 22_500_000, 23_500_000],
           'IL18': ['chr11', 111_500_000, 112_500_000],
           'NLRP1': ['chr17', 5_000_000, 6_000_000],
           'STAT1': ['chr2', 190_500_000, 191_500_000],
           'NFKB1': ['chr4', 102_000_000, 103_000_000],
           'IL1B': ['chr2', 112_250_000, 113_250_000],
           'PTX3': ['chr3', 157_250_000, 158_250_000],
           'CXCL8': ['chr4', 73_250_000, 74_250_000]
          }
"""
regions = {#'group1': ['chr1', 204000000, 207000000],#6h mix
           #'group2': ['chr1', 227000000, 230000000],#6h mix
           'NLRP1': ['chr17', 5_000_000, 6_000_000],
           'group3': ['chr10', 113000000, 115000000]#6h mix
           #'group4': ['chr11', 16000000, 19000000],#6h mix
           #'group5': ['chr4', 73200000, 76200000],#6h mix
           #'group6': ['chr9', 135000000, 138000000],#6h mix
           #'group7': ['chr1', 228000000, 230000000],#24h up
           #'group8': ['chr17', 40000000, 42000000],#24 down
           #'group9': ['chr19', 42700000, 45700000],#24 down
           #'group10': ['chr3', 50000000, 53000000],#24 down
           #'group11': ['chr6', 31300000, 32300000]
          }

figures = []
for gene, region in regions.items():
    norm = LogNorm(vmax=0.1, vmin=0.001)
    data = clr_balanced.matrix(balance=True).fetch(tuple(region))

    fig, ax = plt.subplots(figsize=(20, 10))
    im = pcolormesh_45deg(ax, data, start=region[1], resolution=resolution, norm=norm, cmap='fall')
    #region_length = region[2] - region[1]
    ax.set_aspect(0.5)
    ax.set_ylim(0, 20*windows[0])
    format_ticks(ax, rotate=False)
    ax.xaxis.set_visible(False)

    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="1%", pad=0.1, aspect=6)
    plt.colorbar(im, cax=cax)

    insul_region = bioframe.select(insulation_table, region)
    print(insul_region)
    sys.exit()
    ins_ax = divider.append_axes("bottom", size="50%", pad=0., sharex=ax)

    ins_ax.plot(insul_region[['start', 'end']].mean(axis=1),
                insul_region[f'log2_insulation_score_{windows[0]}'], label=f'Window {windows[0]} bp')

    boundaries = insul_region[~np.isnan(insul_region[f'boundary_strength_{windows[0]}'])]
    weak_boundaries = boundaries[~boundaries[f'is_boundary_{windows[0]}']]
    strong_boundaries = boundaries[boundaries[f'is_boundary_{windows[0]}']]
    ins_ax.scatter(weak_boundaries[['start', 'end']].mean(axis=1),
                weak_boundaries[f'log2_insulation_score_{windows[0]}'], label='Weak boundaries')
    ins_ax.scatter(strong_boundaries[['start', 'end']].mean(axis=1),
                strong_boundaries[f'log2_insulation_score_{windows[0]}'], label='Strong boundaries')

    ins_ax.legend(bbox_to_anchor=(0., -1), loc='lower left', ncol=4)

    format_ticks(ins_ax, y=False, rotate=False)
    ax.set_title(f"{sample} Insulation scores\n{gene}, {region[0]}:{region[1]}-{region[2]}", fontsize=16, pad=20)
    ax.set_xlim(region[1], region[2])
    ins_ax.axhline(y=0, color='black', linestyle='-')
    ins_ax.set_ylim(-1.5, 1.5)

    plt.tight_layout()
    figures.append(fig)

pdf_output = f'{Path(sys.argv[1]).parent.resolve()}/{sample}_insulation_groups_balanced2.pdf'
with PdfPages(pdf_output) as pdf:
    for fig in figures:
        pdf.savefig(fig)
        plt.close(fig)
