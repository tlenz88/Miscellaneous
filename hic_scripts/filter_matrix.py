import numpy
import pandas as pd
import sys
from functools import reduce

mtx = pd.read_csv(sys.argv[1], sep="\t", header=0)
mtx["bin1"] = mtx["bin1"] // 10000 + 1
mtx["bin2"] = mtx["bin2"] // 10000 + 1
print(mtx)

genes = pd.read_csv(sys.argv[2], sep="\t", header=None)
genes[3] = genes[3] // 10000 + 1
genes[4] = genes[4] // 10000 + 1
print(genes)

df = pd.merge(genes, left_on=[0,])
"""
sizes = pd.read_csv(sys.argv[3], sep="\t", header=None)
sizes[1] = sizes[1] // 10000 + 1

start_index = 0
chrom_starts = {}
for chrom in sizes.itertuples():
	chrom_starts[chrom[1]] = start_index
	start_index += chrom[2]

bins = []
for c in chrom_starts:
	chrom_sizes = sizes[sizes[0] == c]
	start = chrom_starts[c]
	for i in genes.itertuples():
		if i[1] == c:
			if i[7] == "+":
				bins.append(i[4] + start)
			elif i[7] == "-":
				bins.append(i[5] + start)

bins = list(dict.fromkeys(bins))

mtx = mtx.loc[mtx["bin1"].isin(bins), : ]
mtx = mtx.loc[mtx["bin2"].isin(bins), : ]
"""
#mtx.to_csv(sys.argv[4], sep="\t", header=False, index=False)