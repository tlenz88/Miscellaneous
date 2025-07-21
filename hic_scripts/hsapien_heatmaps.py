#!/usr/bin/env python3

import sys
import csv
import argparse
import numpy as np
import pandas as pd
import heapq
import math
from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

def parse_args(args):
	parser = argparse.ArgumentParser(description="Check the help flag")
	parser.add_argument("-ch",
						"--chr_sizes",
						dest="sizes",
						help="REQUIRED: chromosome sizes (.chrom.sizes)",
						required=True)
	parser.add_argument("-m",
						"--matrix",
						dest="matrix",
						help="REQUIRED: tab-delimited matrix output from HiC-Pro.",
						required=True)
	parser.add_argument("-s",
						"--sample",
						dest="sample",
						help="Sample name",
						required=False)
	parser.add_argument("-o",
						"--output",
						dest="output",
						help="Output directory.",
						required=False)
	parser.add_argument("-r",
						"--resolution",
						dest="res",
						help="Binning resolution.",
						required=False,
						default=10000)
	parser.add_argument("-n",
						"--normalize",
						dest="norm",
						help="Perform counts-per-million normalization.",
						required=False,
						action="store_false",
						default=None)
	parser.add_argument("-g",
						"--gene_list",
						dest="genes",
						help="List of genes to annotate in gff format.",
						required=False)
	parser.add_argument("-c",
						"--centromeres",
						dest="centromeres",
						help="List of centromere locations.",
						required=False)
	return parser.parse_args()

def input_params(args):
	if args.output is None:
		out = os.path.dirname(args.input)
	else:
		out = args.output
	if args.sample is None:
		sample = os.path.split(os.path.dirname(os.path.abspath(args.matrix)))[-1]
	else:
		sample = args.sample
	if args.res is not None:
		res = resolution(args.res)
	return out, sample, res

def resolution(res):
	res = res.lower()
	if "kb" in res:
		return int(res.split("kb")[0]) * 1000
	elif "mb" in res:
		return int(res.split("mb")[0]) * 1000000
	else:
		return int(res)

def get_chrom_starts(sizes):
	start_index = 0
	chrom_starts = {}
	for chrom in sizes.itertuples():
		chrom_starts[chrom[1]] = start_index
		start_index += chrom[2]
	return chrom_starts

def main():
	args = parse_args(sys.argv[1:])
	out, sample, res = input_params(args)
#	maxColor = []
	sizes = pd.read_csv(args.sizes, sep='\t', header=None)
	sizes[1] = sizes[1] // res + 1
#
#	chrom_starts = get_chrom_starts(sizes)
#
#	df = pd.read_csv(args.matrix, sep='\t', header=None)
#	read_sum = np.sum(df[2])
#	df[2] = df[2] / read_sum * 1000000
#
#	chrom_num = 1
#	for c in chrom_starts:
#		chrom_sizes = sizes[sizes[0] == c]
#		start = chrom_starts[c]
#		end = start + int(chrom_sizes[1])
#
#		sdf = df[df[0] > start]
#		sdf = sdf[sdf[0] < end]
#		sdf = sdf[sdf[1] > start]
#		sdf = sdf[sdf[1] < end]
#		sdf[0] = sdf[0] - start
#		sdf[1] = sdf[1] - start
#
#		arr_size = end - start
#
#		top10percent = int(math.ceil((arr_size**2) * 0.10))
#		mC = min(heapq.nlargest(top10percent, sdf[2].tolist()))
#
#		maxColor.append(mC)
#
#		arr = np.zeros((arr_size,arr_size))
#		for i in sdf.itertuples():
#			if i[1] == i[2]:
#				pass
#			else:
#				arr[int(i[1]) - 1][int(i[2]) - 1] = i[3]
#				arr[int(i[2]) - 1][int(i[1]) - 1] = i[3]
#		np.save("/mnt/f/toxo_project/matrices/NI_HFF/%s.npy" % c, arr)
#
#	del df, sdf, arr, start, end, chrom_sizes, sizes, res, chrom_starts

	sizes = sizes[0].tolist()
	chrom_num = 1
	maxColor = [0.0014951233714438853, 0.002392201896925347, 0.0026479429319667363, 0.0022925740328093376, 0.0016692770203245897, 0.002536458181340988, 0.002699760527757716, 0.0026300000106725683, 0.0026997155016064128, 0.0012164490166722667, 0.0024240428902549817, 0.002778698875350086, 0.0025954499439062146, 0.0025472719620122115, 0.0029267673739093642, 0.00030986997326591306, 0.0017689424062333516, 0.0028219765111106333, 0.0034678766515493895, 0.0027330798797217926, 0.003120087154527963, 0.0036085233395022555, 0.00585935062569534, 0.005634527547881879]
	pdf = PdfPages('%s/%s_10kb.pdf' % (out, sample))
	for c in sizes:
		print(c)
		arr = np.load("/mnt/f/toxo_project/matrices/NI_HFF/%s.npy" % c)
		arr = arr.astype("float32")
		plt.imshow(arr, cmap='Reds', origin='lower', vmax=maxColor[chrom_num-1])
		cbar = plt.colorbar(cmap = 'Reds')
		cbar.ax.set_ylabel('Normalized contact counts', rotation=270, labelpad=12, size=10)
		plt.title("%s chr%s %skb \n intrachromosomal interactions" % (sample, chrom_num, res // 1000), pad=12)
		chrom_num += 1
		pdf.savefig()
		plt.close()
		del arr, cbar
	pdf.close()

if __name__ == "__main__":
	main()