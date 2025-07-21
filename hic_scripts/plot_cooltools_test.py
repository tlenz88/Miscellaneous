import sys
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os, subprocess
import cooler
import cooltools.lib.plotting
import cooltools
import bioframe
from matplotlib.colors import LogNorm
from mpl_toolkits.axes_grid1 import make_axes_locatable

clr = cooler.Cooler(sys.argv[1])
bins = clr.bins()[:]
hg38_genome = bioframe.load_fasta('/mnt/f/organism_genome/Hsapien_autosomes/GRCh38.fasta')
gc_cov = bioframe.frac_gc(bins[['chrom', 'start', 'end']], hg38_genome)

view_df = pd.DataFrame({'chrom': clr.chromnames,
                        'start': 0,
                        'end': clr.chromsizes.values,
                        'name': clr.chromnames}
                      )

cis_eigs = cooltools.eigs_cis( 
                        clr, 
                        gc_cov, 
                        view_df=view_df, 
                        n_eigs=3,
                        )

eigenvector_track = cis_eigs[1][['chrom','start','end','E1']]

f, ax = plt.subplots(
    figsize=(15, 10),
)

norm = LogNorm(vmax=0.1)

im = ax.matshow(
    clr.matrix()[:], 
    norm=norm,  
    cmap='fall'
); 
plt.axis([0,500,500,0])

divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.1)
plt.colorbar(im, cax=cax, label='corrected frequencies');
ax.set_ylabel('chr2:0-50Mb')
ax.xaxis.set_visible(False)

ax1 = divider.append_axes("top", size="20%", pad=0.25, sharex=ax)
weights = clr.bins()[:]['weight'].values
ax1.plot([0,500],[0,0],'k',lw=0.25)
ax1.plot( eigenvector_track['E1'].values, label='E1')

ax1.set_ylabel('E1')
ax1.set_xticks([])


for i in np.where(np.diff( (cis_eigs[1]['E1']>0).astype(int)))[0]:
    ax.plot([0, 500],[i,i],'k',lw=0.5)
    ax.plot([i,i],[0, 500],'k',lw=0.5)

plt.savefig(f'compartment_plot_{chrom}.pdf', bbox_inches='tight')
plt.close()

