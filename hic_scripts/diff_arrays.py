#!/usr/bin/env python3

import sys
import csv
import argparse
import numpy as np
import pandas as pd
import heapq
import re
import math
import xarray as xr
from matplotlib import pyplot as plt
from matplotlib import colors as colors
from matplotlib.backends.backend_pdf import PdfPages

def parse_args(args):
	parser = argparse.ArgumentParser(description="Check the help flag")
	parser.add_argument("-s",
						"--sample",
						dest="sample",
						help="REQUIRED: sample names",
						required=False)
	parser.add_argument("-c",
						"--chr_sizes",
						dest="sizes",
						help="REQUIRED: chromosome sizes (.chrom.sizes)",
						required=True)
	parser.add_argument("-m1",
						"--matrix1",
						dest="matrix1",
						help="REQUIRED: tab-delimited contact count matrix for experimental sample",
						required=True)
	parser.add_argument("-m2",
						"--matrix2",
						dest="matrix2",
						help="REQUIRED: tab-delimited contact count matrix for control sample",
						required=True)
	parser.add_argument("-o",
						"--output",
						dest="out",
						help="REQUIRED: output directory",
						required=True)
	parser.add_argument("-r",
						"--resolution",
						dest="res",
						help="REQUIRED: resolution",
						required=True)
	return parser.parse_args()

def resolution(res):
	res = res.lower()
	if res.isnumeric():
		return res
	elif "kb" in res:
		return int(res.split("kb")[0]) * 1000
	elif "mb" in res:
		return int(res.split("mb")[0]) * 1000000

def chrom(bin_counts, res):
	bins = {}
	chr = []
	num = 0
	c = 0
	for i in bin_counts:
		num += int(i[1])
		chr.append(i[0])
		bins[i[0]] = {}
		for j in range(int(i[1])):
			mid = int(res / 2) + j * res
			bins[i[0]][mid] = c
			c += 1
	return bins, chr, num, c

def read_chrom_sizes(bin_counts):
	chrom_sizes = {}
	for i in bin_counts:
		chrom_sizes[i[0]] = int(i[1])
	return chrom_sizes

def get_chrom_starts(chrom_sizes, chr):
	start_index = 0
	end_index = 0
	chrom_starts = {}
	chrom_ends = {}
	for chrom in chr:
		end_index = start_index + chrom_sizes[chrom]
		chrom_starts[chrom] = start_index
		chrom_ends[chrom] = end_index
		start_index = end_index
	return chrom_starts

def plot_diff_intra(pdf,m,maxColor,sample,ch,res):
	norm = colors.TwoSlopeNorm(vmin=-maxColor, vcenter=0, vmax=maxColor)
	plt.imshow(m, cmap = 'RdBu_r', origin = 'lower', norm=norm)
	cbar = plt.colorbar(cmap = 'RdBu_r')
	cbar.ax.set_ylabel('log$_\mathrm{2}$ fold change', rotation = 270, labelpad = 12, size = 10)
	plt.title("%s chr%s %skb\ndifferential intrachromosomal interactions" % (sample, ch, res // 1000))
	pdf.savefig()
	plt.close()

def plot_diff_inter(pdf,m_inter,maxColor,sample,res):
	norm = colors.TwoSlopeNorm(vmin=-maxColor, vcenter=0, vmax=maxColor)
	plt.imshow(m_inter, cmap = 'RdBu_r', origin = 'lower', norm=norm)
	cbar = plt.colorbar(cmap = 'RdBu_r')
	cbar.ax.set_ylabel('log$_\mathrm{2}$ fold change', rotation = 270, labelpad = 12, size = 10)
	plt.title("%s %skb\ninterchromosomal interactions" % (sample, res // 1000))
	pdf.savefig()
	plt.close()

def main():
	args = parse_args(sys.argv[1:])
	res = int(resolution(args.res))
	with open(args.sizes, 'r') as f:
		bin_counts = []
		for i in csv.reader(f, delimiter='\t'):
			i[1] = int(i[1]) // res + 1
			bin_counts.append(i)
	bins, chr, num, c = chrom(bin_counts, res)
	chrom_sizes = read_chrom_sizes(bin_counts)
	chrom_starts = get_chrom_starts(chrom_sizes, chr)

	mat1 = np.zeros((num, num))
	with open(args.matrix1, 'r') as mtx:
		for i in csv.reader(mtx, delimiter = "\t"):
			mat1[int(i[0]) - 1][int(i[1]) - 1] = i[2]
			mat1[int(i[1]) - 1][int(i[0]) - 1] = i[2]
		reads1 = np.sum(mat1)
		mat1 = mat1 / (reads1 * 0.000001)
	for i in range(0,len(mat1)-2):
		mat1[i,i] = 0

	mat2 = np.zeros((num, num))
	with open(args.matrix2, 'r') as mtx:
		for i in csv.reader(mtx, delimiter = "\t"):
			mat2[int(i[0]) - 1][int(i[1]) - 1] = i[2]
			mat2[int(i[1]) - 1][int(i[0]) - 1] = i[2]
		reads2 = np.sum(mat2)
		mat2 = mat2 / (reads2 * 0.000001)
	for i in range(0,len(mat2)-2):
		mat2[i,i] = 0

	if args.sample:
		sample = str(args.sample)
		pdf = PdfPages('%s/%s_%skb_diff_int.pdf' % (str(args.out), sample, res // 1000))
	else:
		sample = ""
		pdf = PdfPages('%s/%skb_diff_int.pdf' % (str(args.out), res // 1000))
	max_vals = []
	maxVals = [1.2760485514273467, 1.2095667276677873, 1.42545629878426, 1.3688112667713461, 
	1.550765399891101, 1.4799710083212327, 1.635097195451539, 1.5287204436767086, 
	1.5430962423708658, 1.5596924300348314, 1.584426206190232, 1.6011398124472545, 
	1.5935503177027448, 1.615250104158752]
	chrom_color = 0
	for c in chrom_starts:
		start = chrom_starts[c]
		end = start + chrom_sizes[c]
		m1 = mat1[start:end, start:end]
		m2 = mat2[start:end, start:end]
		m = np.log2(np.true_divide(m1, m2, out=np.zeros_like(m1), where=m2!=0))
		m = np.nan_to_num(m, neginf=0)
		allContacts = []
		for i, row in enumerate(m):
			for j, elem in enumerate(row):
				allContacts.append(abs(elem))
		top10percent = int(math.ceil(len(allContacts) * 0.10))
		#maxColor = min(heapq.nlargest(top10percent, allContacts))
		maxColor = maxVals[chrom_color]
		chrom_color += 1
		#m[np.abs(m) < 1] = 0
		if len(str(c)) > 2:
			ch = re.split('_', str(c))[1]
		else:
			ch = str(c)
		plot_diff_intra(pdf,m,maxColor,sample,ch,res)
	print(max_vals)
	for key, val in chrom_starts.items():
		mat1[chrom_starts[key]:(chrom_starts[key]+chrom_sizes[key]), chrom_starts[key]:(chrom_starts[key]+chrom_sizes[key])] = 0
		mat2[chrom_starts[key]:(chrom_starts[key]+chrom_sizes[key]), chrom_starts[key]:(chrom_starts[key]+chrom_sizes[key])] = 0
	m_inter = np.log2(np.true_divide(mat2, mat1, out=np.zeros_like(mat2), where=mat1!=0))
	m_inter = np.nan_to_num(m_inter, neginf=0)
	allContacts = []
	for i, row in enumerate(m_inter):
		for j, elem in enumerate(row):
			allContacts.append(abs(elem))
	top10percent = int(math.ceil(len(allContacts) * 0.10))
	maxColor = min(heapq.nlargest(top10percent, allContacts))
	plot_diff_inter(pdf,m_inter,maxColor,sample,res)
	pdf.close()

if __name__ == "__main__":
	main()