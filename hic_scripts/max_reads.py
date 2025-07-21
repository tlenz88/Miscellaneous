import sys
import csv
import argparse
import numpy as np
import pandas as pd
import heapq
import re
import math

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

def main():
	res = 10000
	with open(sys.argv[1], 'r') as f:
		bin_counts = []
		for j in csv.reader(f, delimiter='\t'):
			j[1] = int(j[1]) // res + 1
			bin_counts.append(j)
	bins, chr, num, c = chrom(bin_counts, res)
	chrom_sizes = read_chrom_sizes(bin_counts)
	chrom_starts = get_chrom_starts(chrom_sizes, chr)
	mat = np.zeros((num, num))
	df = np.zeros((14,len(sys.argv)-2))
	gmax = []
	sample = 0
	for i in range(2,len(sys.argv)):
		chromosome = 0
		with open(sys.argv[i], 'r') as mtx:
			for k in csv.reader(mtx, delimiter = "\t"):
				if int(k[0]) <= 2337 and int(k[1]) <= 2337:
					mat[int(k[0]) - 1][int(k[1]) - 1] = k[2]
					mat[int(k[1]) - 1][int(k[0]) - 1] = k[2]
			reads = np.sum(mat)
		for j in range(0,len(mat)-2):
			mat[j,j] = 0
		for c in chrom_starts:
			start = chrom_starts[c]
			end = start + chrom_sizes[c]
			m = mat[start:end, start:end]
			allContacts = []
			for k, row in enumerate(m):
				for l, elem in enumerate(row):
					allContacts.append(elem)
			top10percent = int(math.ceil(len(allContacts) * 0.10))
			maxColor = min(heapq.nlargest(top10percent, allContacts))
			df[chromosome,sample] = maxColor
			chromosome += 1
		for key, val in chrom_starts.items():
			mat[chrom_starts[key]:(chrom_starts[key]+chrom_sizes[key]),
			chrom_starts[key]:(chrom_starts[key]+chrom_sizes[key])] = 0
		allContacts = []
		for j, row in enumerate(mat):
			for k, elem in enumerate(row):
				allContacts.append(elem)
		top10percent = int(math.ceil(len(allContacts) * 0.10))
		maxColor = min(heapq.nlargest(top10percent, allContacts))
		gmax.append(maxColor)
		sample += 1
	print(df.max(axis=1).tolist())
	print(max(gmax))

if __name__ == "__main__":
	main()