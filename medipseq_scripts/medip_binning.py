import sys
import csv
import re
import pandas as pd
import numpy as np
import os

def find_genes(gff):
	# split gff df by type
	genes = gff[(gff[2] == "gene") | (gff[2] == "ncRNA_gene")].copy()

	# extract accession numbers for genes and exons
	genes[9] = genes[8].str.extract(r'ID=(.*?);', expand=True)

	# drop unnecessary columns from dfs
	genes.drop(columns=[1,2,5,7,8], inplace=True)

	# sort columns to fit .bed format (for downstream applications)
	genes = genes[[0,3,4,9,10,6]]

	# set column names
	genes.columns = ["chr", "start", "end", "gene", "score", "strand"]

	return genes


def find_exons(gff, genes):
	exons = gff[gff[2] == "exon"].copy()

	# extract accession numbers for genes and exons
	exons[9] = exons[8].str.extract(r'ID=exon_(.*?)\.', expand=True)
	exons[11] = exons[8].str.extract(r'\.\d-(.*?);', expand=True)
	exons[12] = exons[8].str.extract(r'\.(\d)-', expand=True)
	exons[13] = exons[8].str.extract(r'\.(\d);', expand=True)

	# filter list by 
	exons[12] = exons[12].astype(int)
	exons[13] = exons[13].astype(int)
	exons = exons[(exons[12] == 1) | (exons[13] == 1)]
	

	# drop unnecessary columns from df
	exons.drop(columns=[1,2,5,7,8,12,13], inplace=True)

	# sort columns to fit .bed format (for downstream applications)
	exons = exons[[0,3,4,9,10,6,11]]

	# set column names
	exons.columns = ["chr", "start", "end", "gene", "score", "strand", "exon_num"]

	# filter exons not in protein coding genes
	gene_list = list(genes["gene"])
	exons = exons[exons["gene"].isin(gene_list)]
	exons.reset_index(inplace = True, drop=True)

	return exons


def find_introns(exons):
	# create empty df to store intron coordinates
	introns = pd.DataFrame(columns=["chr", "start", "end", "gene", "score", "strand", "exon_num"])
	current_gene = ""
	exon_num = ""
	start = 0

	# iterate over exons to find introns
	for i in exons.itertuples():
		if i[4] == current_gene:
			if i[6] == "-":
				exon_num = i[7]
			current_intron = {"chr":i[1], "start":start, "end":i[2]-1, "gene":i[4], "score":1, "strand":i[6], "exon_num":exon_num}
			introns = pd.concat([introns, pd.DataFrame([current_intron])], ignore_index=True)
			start = i[3] + 1
			exon_num = i[7]
		else:
			current_gene = i[4]
			start = i[3] + 1
			exon_num = i[7]

	return introns


def find_five_prime(genes):
	five_prime = genes.copy()
	starts = []
	ends = []
	for i in five_prime.itertuples():
		if i[6] == "+":
			starts.append(i[2] - 1000)
			ends.append(i[2] - 1)
		else:
			starts.append(i[3] + 1)
			ends.append(i[3] + 1000)
	five_prime["start"] = starts
	five_prime["end"] = ends
	return five_prime


def find_three_prime(genes):
	three_prime = genes.copy()
	starts = []
	ends = []
	for i in three_prime.itertuples():
		if i[6] == "+":
			starts.append(i[3] + 1)
			ends.append(i[3] + 1000)
		else:
			starts.append(i[2] - 1000)
			ends.append(i[2] - 1)
	three_prime["start"] = starts
	three_prime["end"] = ends
	return three_prime


def binning(df, feature):
	binned_feature = []
	for i in feature.chr.unique():
		sub_df = df.loc[df.chr == i]
		sub_feature = feature.loc[feature.chr == i]
		for j in sub_feature.itertuples():
			feature_bins = sub_df.loc[(sub_df.pos >= j[2]) & (sub_df.pos <= j[3])].copy()
			feature_bins = feature_bins.reset_index(drop=True)
			binsize = len(feature_bins) / 5
			if j[6] == "+":
				bin1 = feature_bins["reads"].iloc[:round(binsize)].sum()
				bin2 = feature_bins["reads"].iloc[round(binsize):round(binsize*2)].sum()
				bin3 = feature_bins["reads"].iloc[round(binsize*2):round(binsize*3)].sum()
				bin4 = feature_bins["reads"].iloc[round(binsize*3):round(binsize*4)].sum()
				bin5 = feature_bins["reads"].iloc[round(binsize*4):].sum()
			else:
				bin5 = feature_bins["reads"].iloc[:round(binsize)].sum()
				bin4 = feature_bins["reads"].iloc[round(binsize):round(binsize*2)].sum()
				bin3 = feature_bins["reads"].iloc[round(binsize*2):round(binsize*3)].sum()
				bin2 = feature_bins["reads"].iloc[round(binsize*3):round(binsize*4)].sum()
				bin1 = feature_bins["reads"].iloc[round(binsize*4):].sum()
			try:
				binned_feature.append([j[4], j[7], binsize, bin1, bin2, bin3, bin4, bin5])
			except:
				binned_feature.append([j[4], binsize, bin1, bin2, bin3, bin4, bin5])
	return pd.DataFrame(binned_feature)


def zscore_by_row(df):
	df[8] = df[[3,4,5,6,7]].mean(axis=1)
	df[9] = df[[3,4,5,6,7]].std(axis=1)
	for i in range(3,7):
		df[i+7] = (df[i] - df[8]) / df[9]
		df.loc[(df[i] == 0), i+7] = "NA"


def main():
	# read input bed file
	df = pd.read_csv(sys.argv[1], sep = '\t', header = None, names = ("chr", "pos", "reads"))

	# generate output file paths/names
	outdir = os.path.dirname(os.path.abspath(sys.argv[1]))
	genes_output = "".join([os.path.splitext(os.path.basename(sys.argv[1]))[0], "_genes_binned.txt"])
	exons_output = "".join([os.path.splitext(os.path.basename(sys.argv[1]))[0], "_exons_binned.txt"])
	introns_output = "".join([os.path.splitext(os.path.basename(sys.argv[1]))[0], "_introns_binned.txt"])
	five_prime_output = "".join([os.path.splitext(os.path.basename(sys.argv[1]))[0], "_five_prime_binned.txt"])
	three_prime_output = "".join([os.path.splitext(os.path.basename(sys.argv[1]))[0], "_three_prime_binned.txt"])

	# read gff file
	gff = pd.read_csv(sys.argv[2], sep = '\t', header = None)

	# add column to act as score for .bed reformat
	gff[10] = 1

	# obtain features from gff
	genes = find_genes(gff)
	exons = find_exons(gff, genes)
	introns = find_introns(exons)
	five_prime = find_five_prime(genes)
	three_prime = find_three_prime(genes)

	# bin df reads within features
	binned_genes = binning(df, genes)
	binned_exons = binning(df, exons)
	binned_introns = binning(df, introns)
	binned_five_prime = binning(df, five_prime)
	binned_three_prime = binning(df, three_prime)

	# write binned feature counts to tab-delimited text file
	binned_genes.to_csv("".join([outdir, "/", genes_output]), sep="\t", index=False, header=False)
	binned_exons.to_csv("".join([outdir, "/", exons_output]), sep="\t", index=False, header=False)
	binned_introns.to_csv("".join([outdir, "/", introns_output]), sep="\t", index=False, header=False)
	binned_five_prime.to_csv("".join([outdir, "/", five_prime_output]), sep="\t", index=False, header=False)
	binned_three_prime.to_csv("".join([outdir, "/", three_prime_output]), sep="\t", index=False, header=False)

if __name__ == "__main__":
	main()