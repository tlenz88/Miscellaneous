#!/usr/bin/env python

import sys
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import colors as colors
from matplotlib import font_manager
import matplotlib.font_manager as fm

genes = pd.read_csv(sys.argv[2], sep='\t', header=None)
genes = genes[genes[2] == 'protein_coding_gene']
genes[3] = round(genes[3], -4)
genes[4] = round(genes[4], -4)
gene_bins1 = pd.DataFrame({
    'chr1': genes[0],
    'bin1': genes.apply(lambda x: x[3] if x[6] == '+' else x[4], axis=1),
    'gene1': range(1, len(genes) + 1)
})
gene_bins2 = pd.DataFrame({
    'chr2': genes[0],
    'bin2': genes.apply(lambda x: x[3] if x[6] == '+' else x[4], axis=1),
    'gene2': range(1, len(genes) + 1)
})

diff = pd.read_csv(sys.argv[1], sep='\t', names=["chr1", "bin1", "chr2", "bin2", "pval", "log2fc"])
diff_gene = pd.merge(diff, gene_bins1, on=['chr1', 'bin1'])
diff_gene = pd.merge(diff_gene, gene_bins2, on=['chr2', 'bin2'])

diff_arr = np.zeros((len(genes), len(genes)))
for i in diff_gene.itertuples():
    diff_arr[int(i[7]) - 1][int(i[8]) - 1] = i[6]
    diff_arr[int(i[8]) - 1][int(i[7]) - 1] = i[6]

chr_sizes = pd.read_csv(sys.argv[3], sep='\t')
chr_colors = ['#00FF00', '#FF00FF', '#007FFF', '#FF7F00', '#7FBF7F', 
              '#6D258F', '#FF0000', '#038021', '#16FCDF', '#F7FE1E', 
              '#0000FF', '#FF7FFF', '#FD4280', '#904E0E']

maxColor = max(np.max(diff_arr), abs(np.min(diff_arr)))

chr_counts = pd.DataFrame(genes[0].value_counts()).sort_values(0).reset_index()
chr_counts['start'] = 0
for i in range(len(chr_counts)):
    if i == 0:
        chr_counts.at[i, 'start'] = chr_counts.at[i, 'count']
    else:
        chr_counts.at[i, 'start'] = chr_counts.at[i - 1, 'start'] + chr_counts.at[i, 'count']
chr_counts['start'] = chr_counts['start'].shift(periods=1, fill_value=0)

pdf = PdfPages(sys.argv[4])
fig, ax = plt.subplots()
norm = colors.TwoSlopeNorm(vmin=-maxColor, vcenter=0, vmax=maxColor)
plt.imshow(diff_arr, cmap='RdBu_r', origin='lower', norm=norm)
#fontprop = fm.FontProperties(fname='/usr/share/fonts/truetype/msttcorefonts/Times_New_Roman.ttf')
plt.rcParams['font.size'] = 8
ax.set_xticks([])
ax.set_yticks([])
plt.gca().set_aspect('equal')
for axis in ['top', 'bottom', 'left', 'right']:
    ax.spines[axis].set_linewidth(1)
cbar = plt.colorbar(cmap='RdBu_r')
cbar.ax.set_ylabel(r'log$_\mathrm{2}$ fold change', rotation=270, labelpad=10, size=8) #, fontproperties=fontprop)
cbar.ax.tick_params(labelsize=6)
ax.set(xlim=(-0.5, 60.5), ylim=(-0.5, 60.5))
for c in chr_counts.itertuples():
    ax.arrow(c[3] - 0.4, -1, c[2], 0, 
    width=1, head_width=0, head_length=0, 
    facecolor=chr_colors[c[0]], edgecolor='none', 
    length_includes_head=True, clip_on=False)
    plt.text(c[3] + 0.5*c[2] + 0.5, -2, c[1], ha='right', va='top', rotation=45, size=6) #, fontproperties=fontprop)
    ax.arrow(-1, c[3] - 0.4, 0, c[2], 
    width=1, head_width=0, head_length=0, 
    facecolor=chr_colors[c[0]], edgecolor='none', 
    length_includes_head=True, clip_on=False)
    plt.text(-2, c[3] + 0.5*c[2] - 0.5, c[1], ha='right', va='center', size=6)

plt.tight_layout()
pdf.savefig()
plt.close()
pdf.close()
