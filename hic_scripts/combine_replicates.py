#!/usr/bin/env python3

"""
Created: February 10, 2022
Updated: March 14, 2023
Author(s): Todd Lenz, tlenz001@ucr.edu

Combines HiC-Pro replicate matrices.
"""

import sys
import pandas as pd
import numpy as np
import os
from functools import partial, reduce

def merge_hic():
	df = {}
	for i in range(1,len(sys.argv)):
		if os.path.splitext(sys.argv[i])[1] == ".matrix":
			sample = os.path.dirname(sys.argv[i])
			df[sample] = pd.read_csv(sys.argv[i], sep='\t', names=("b1", "b2", "val%s" % i))
		elif sys.argv[i].isdigit():
			res = sys.argv[i]
	df_final = reduce(partial(pd.merge, on=['b1', 'b2'], how='outer'), df.values()).sort_values(['b1','b2']).fillna(0)
	col_sums = list(df_final.sum()[2:len(df_final.columns)])
	total_reads = sum(col_sums)
	weighted_colnames = []
	for i in range(1,len(df_final.columns)-1):
		df_final["".join(["weighted", str(i)])] = df_final["val%s" % i] * col_sums[i-1] / total_reads
		weighted_colnames.append("".join(["weighted", str(i)]))
	df_final["weighted_sum"] = df_final[weighted_colnames].sum(axis=1)
	df_final = df_final[['b1', 'b2', 'weighted_sum']]
	try:
		outfile = "".join([os.path.split(os.path.split(os.path.abspath(sys.argv[1]))[0])[0], "/samples_merged_", res, ".matrix"])
	except:
		outfile = "".join([os.path.split(os.path.split(os.path.abspath(sys.argv[1]))[0])[0], "/samples_merged.matrix"])
	df_final.to_csv(outfile, sep = '\t', index = False, header = False)

def main():
	df = {}
	for i in range(1,len(sys.argv)):
		if os.path.splitext(sys.argv[i])[1] == ".bed":
			sample = os.path.dirname(sys.argv[i])
			df[sample] = pd.read_csv(sys.argv[i], sep='\t', names=("chr", "coord", "val%s" % i))
		elif sys.argv[i].isdigit():
			res = sys.argv[i]
	df_final = reduce(partial(pd.merge, on=["chr", "coord"], how='outer'), df.values()).sort_values(["chr", "coord"]).fillna(0)
	col_sums = list(df_final.sum()[2:len(df_final.columns)])
	total_reads = sum(col_sums)
	weighted_colnames = []
	for i in range(1,len(df_final.columns)-1):
		df_final["".join(["weighted", str(i)])] = df_final["val%s" % i] * col_sums[i-1] / total_reads
		weighted_colnames.append("".join(["weighted", str(i)]))
	df_final["weighted_sum"] = df_final[weighted_colnames].sum(axis=1)
	df_final = df_final[["chr", "coord", "weighted_sum"]]
	try:
		outfile = "".join([os.path.split(os.path.split(os.path.abspath(sys.argv[1]))[0])[0], "/samples_merged_", res, ".matrix"])
	except:
		outfile = "".join([os.path.split(os.path.split(os.path.abspath(sys.argv[1]))[0])[0], "/samples_merged.matrix"])
	df_final.to_csv(outfile, sep = '\t', index = False, header = False)

if __name__ == "__main__":
	main()