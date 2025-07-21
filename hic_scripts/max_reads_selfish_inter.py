import sys
import csv
import numpy as np
import pandas as pd
import heapq
import math
import xarray as xr

res = 10000
sizes = pd.read_csv(sys.argv[1], sep='\t', header=None)
sizes[2] = sizes[1] // res + 1

max_vals = []
for df in sys.argv[2:]:
	df = pd.read_csv(df, sep='\t')
	df["bin1"] = df["bin1"] // res
	df["bin2"] = df["bin2"] // res

	df_inter = df[df.chr1 != df.chr2]
	if len(df_inter) > 1:
		sizes_inter = sizes.copy()
		starts = []
		ends = []
		start_index = 0
		for i in sizes_inter.itertuples():
			starts.append(start_index)
			start_index += i[3]
		sizes_inter[1] = starts
		sizes_inter[2] = sizes_inter[2] + sizes_inter[1] - 1

		sizes_inter.columns = ["chr1", "start", "end"]
		df_inter = df[df.chr1 != df.chr2]
		df_inter = df_inter.merge(sizes_inter, on="chr1")
		df_inter["bin1"] = df_inter["bin1"] + df_inter["start"]
		df_inter.drop(["start", "end"], axis=1, inplace=True)

		sizes_inter.columns = ["chr2", "start", "end"]
		df_inter = df_inter.merge(sizes_inter, on="chr2")
		df_inter["bin2"] = df_inter["bin2"] + df_inter["start"]
		df_inter.drop(["start", "end"], axis=1, inplace=True)
		df_inter = df_inter[df_inter.chr1 != df_inter.chr2]

		max_vals.append(min(heapq.nlargest(math.ceil(len(df_inter) * 0.10), abs(df_inter.log2FC))))

print(max(max_vals))
