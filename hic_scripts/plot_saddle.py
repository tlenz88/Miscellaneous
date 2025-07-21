import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import os

saddle = np.load(sys.argv[1], allow_pickle=True)
fig, ax = plt.subplots(figsize=(8, 8))
norm = LogNorm(vmin=np.min(saddle['saddledata']), vmax=np.max(saddle['saddledata']))
im = plt.imshow(saddle['saddledata'], cmap='RdBu_r', norm = norm)
ax.set_xlabel('E1 quantiles')
ax.set_ylabel('E1 quantiles')
plt.colorbar(im, label='obs/exp', pad=0.025, shrink=0.7)
plt.savefig(f'/mnt/f/toxo_project/HiC/Hsapien_output/output_files/{os.path.basename(os.path.splitext(sys.argv[1])[0])}.pdf', bbox_inches='tight')
plt.close()
