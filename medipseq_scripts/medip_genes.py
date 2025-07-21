import sys
import csv
import re
import pandas as pd
import numpy as np
import os
from scipy.stats import zscore

# read input bed file
df = pd.read_csv(sys.argv[1], sep = '\t', header = None, names = ("chr", "pos", "reads"))

# generate output file paths/names
#outdir = os.path.dirname(os.path.abspath(sys.argv[1]))
outdir = sys.argv[3]
binned_output = "".join([os.path.splitext(os.path.basename(sys.argv[1]))[0], "_binned.txt"])
exons_output = "".join([os.path.splitext(os.path.basename(sys.argv[1]))[0], "_exons.txt"])
genes_output = "".join([os.path.splitext(os.path.basename(sys.argv[1]))[0], "_genes.txt"])
intergenic_output = "".join([os.path.splitext(os.path.basename(sys.argv[1]))[0], "_intergenic.txt"])
five_prime_output = "".join([os.path.splitext(os.path.basename(sys.argv[1]))[0], "_five_prime.txt"])
introns_output = "".join([os.path.splitext(os.path.basename(sys.argv[1]))[0], "_introns.txt"])

# read gff files
genes = pd.read_csv(sys.argv[2], sep = '\t', header = None)

# extract exons from gff file
exons = genes[genes[2] == "exon"]
exon_names = exons[8].str.extract(r'(PF3D7_.{7})', expand=True)
exons = exons[[0,3,4]]
exons.columns = ["chr", "start", "end"]
exons["gene"] = exon_names

# extract protein coding genes from gff file
genes = genes[(genes[2] == "protein_coding_gene") | (genes[2] == "ncRNA") | (genes[2] == "rRNA") | (genes[2] == "tRNA") | (genes[2] == "snoRNA") | (genes[2] == "snRNA") | (genes[2] == "pseudogene")]
gene_names = genes[8].str.extract(r'(PF3D7_.{7})', expand=True)
genes = genes[[0,3,4,6]]
genes.columns = ["chr", "start", "end", "strand"]
genes["gene"] = gene_names

# create counters for methylation within features
intergenic_meth = 0

# create empty dataframes to store feature counts
gene_counts = pd.DataFrame(columns=["chr", "pos", "reads", "gene"])
exon_counts = pd.DataFrame(columns=["chr", "pos", "reads", "gene"])
intron_counts = pd.DataFrame(columns=["chr", "pos", "reads", "gene"])
intergenic_counts = pd.DataFrame(columns=["chr", "pos", "reads"])
five_prime_counts = pd.DataFrame(columns=["chr", "pos", "reads", "gene"])
binned_counts = []

# iterate through chromosomes
for i in genes["chr"].unique():

	# create sub dataframes for features within current chromosome
	sub_df = df[df.chr == i]
	sub_genes = genes[genes["chr"] == i]

	# create dataframe from which to filter non-intergenic features
	sub_intergenic = sub_df.copy()

	# iterate through genes within current chromosome
	for j in sub_genes.itertuples():

		# append gene data frame with bins for current gene
		gb = sub_df.loc[(sub_df["pos"] >= j[2]) & (sub_df["pos"] <= j[3])].copy()
		current_gene = str(j[5])
		gene_bins = gb.copy()
		gene_bins["gene"] = current_gene
		gene_counts = pd.concat([gene_counts,gene_bins], ignore_index=True, sort=False)

		if j[4] == "+":
			fp = sub_df.loc[((sub_df["pos"] >= j[2] - 1000) & (sub_df["pos"] < j[2]))].copy()
			gene_5p = fp.copy()
			gene_5p["gene"] = current_gene
			gene_5p = gene_5p.reset_index(drop = True)
			gene_5p1 = gene_5p["reads"].iloc[:200].sum()
			gene_5p2 = gene_5p["reads"].iloc[200:400].sum()
			gene_5p3 = gene_5p["reads"].iloc[400:600].sum()
			gene_5p4 = gene_5p["reads"].iloc[600:800].sum()
			gene_5p5 = gene_5p["reads"].iloc[800:].sum()

			gene_body = gene_bins.reset_index(drop = True)
			binsize = len(gene_body)/5
			gene_body1 = gene_body["reads"].iloc[:round(binsize)].sum()
			gene_body2 = gene_body["reads"].iloc[round(binsize):round(binsize*2)].sum()
			gene_body3 = gene_body["reads"].iloc[round(binsize*2):round(binsize*3)].sum()
			gene_body4 = gene_body["reads"].iloc[round(binsize*3):round(binsize*4)].sum()
			gene_body5 = gene_body["reads"].iloc[round(binsize*4):].sum()

			gene_3p = sub_df.loc[((sub_df["pos"] > j[3]) & (sub_df["pos"] <= j[3] + 1000))]
			gene_3p = gene_3p.reset_index(drop = True)
			gene_3p1 = gene_3p["reads"].iloc[:200].sum()
			gene_3p2 = gene_3p["reads"].iloc[200:400].sum()
			gene_3p3 = gene_3p["reads"].iloc[400:600].sum()
			gene_3p4 = gene_3p["reads"].iloc[600:800].sum()
			gene_3p5 = gene_3p["reads"].iloc[800:].sum()

		else:
			gene_3p = sub_df.loc[((sub_df["pos"] >= j[2] - 1000) & (sub_df["pos"] < j[2]))].copy()
			gene_3p = gene_3p.reset_index(drop = True)
			gene_3p5 = gene_3p["reads"].iloc[:200].sum()
			gene_3p4 = gene_3p["reads"].iloc[200:400].sum()
			gene_3p3 = gene_3p["reads"].iloc[400:600].sum()
			gene_3p2 = gene_3p["reads"].iloc[600:800].sum()
			gene_3p1 = gene_3p["reads"].iloc[800:].sum()

			gene_body = gene_bins.reset_index(drop = True)
			binsize = len(gene_body)/5
			gene_body5 = gene_body["reads"].iloc[:round(binsize)].sum()
			gene_body4 = gene_body["reads"].iloc[round(binsize):round(binsize*2)].sum()
			gene_body3 = gene_body["reads"].iloc[round(binsize*2):round(binsize*3)].sum()
			gene_body2 = gene_body["reads"].iloc[round(binsize*3):round(binsize*4)].sum()
			gene_body1 = gene_body["reads"].iloc[round(binsize*4):].sum()

			fp = sub_df.loc[((sub_df["pos"] > j[3]) & (sub_df["pos"] <= j[3] + 1000))].copy()
			gene_5p = fp.copy()
			gene_5p["gene"] = current_gene
			gene_5p = gene_5p.reset_index(drop = True)
			gene_5p5 = gene_5p["reads"].iloc[:200].sum()
			gene_5p4 = gene_5p["reads"].iloc[200:400].sum()
			gene_5p3 = gene_5p["reads"].iloc[400:600].sum()
			gene_5p2 = gene_5p["reads"].iloc[600:800].sum()
			gene_5p1 = gene_5p["reads"].iloc[800:].sum()

		all_bins = [current_gene, binsize, gene_5p1, gene_5p2, gene_5p3, gene_5p4, gene_5p5, gene_body1, gene_body2, gene_body3, gene_body4, gene_body5, gene_3p1, gene_3p2, gene_3p3, gene_3p4, gene_3p5]
		binned_counts.append(all_bins)

		# drop rows from df containing counts within genes (to get intergenic counts)
		sub_intergenic = sub_intergenic[(sub_intergenic["pos"] < j[2]) | (sub_intergenic["pos"] > j[3])]

		# create sub dataframe for exon(s) within current gene
		sub_exons = exons[exons["gene"] == j[5]]

		current_intron = ""
		intron_start = 0
		intron_end = 0

		# iterate through exons within current gene
		for k in sub_exons.itertuples():

			# create data frame for read counts within current exon
			eb = sub_df.loc[((sub_df["pos"] >= k[2]) & (sub_df["pos"] <= k[3]))]
			exon_bins = eb.copy()
			exon_bins["gene"] = current_gene

			# append methylation counts for current exon
			exon_counts = pd.concat([exon_counts,exon_bins], ignore_index=True, sort=False)

			if current_intron == k[4]:
				intron_end = k[2]
				ib = sub_df.loc[((sub_df["pos"] > intron_start) & (sub_df["pos"] < intron_end))]
				intron_bins = ib.copy()
				intron_bins["gene"] = current_gene
				intron_counts = pd.concat([intron_counts,intron_bins], ignore_index=True, sort=False)
			else:
				current_intron = k[4]
				intron_start = k[3]
			intron_start = k[3]

		# get intergenic read counts
		if int(j[0]) == int(sub_genes.iloc[-1].name):
			intergenic_meth += sub_intergenic["reads"].sum()

		five_prime_counts = pd.concat([five_prime_counts, gene_5p], ignore_index=True, sort=False)

	intergenic_counts = pd.concat([intergenic_counts, sub_intergenic], ignore_index=True, sort=False)

# convert list of lists with 5', gene, and 3' bin counts to a dataframe
binned_counts = pd.DataFrame(binned_counts, columns = ["gene", "bs", "fp1", "fp2", "fp3", "fp4", "fp5", "gb1", "gb2", "gb3", "gb4", "gb5", "tp1", "tp2", "tp3", "tp4", "tp5"])

# drop duplicates from feature count dataframes
gene_counts.drop_duplicates(subset=["chr","pos"], keep="first", inplace=True)
exon_counts.drop_duplicates(subset=["chr","pos"], keep="first", inplace=True)
intron_counts.drop_duplicates(subset=["chr","pos"], keep="first", inplace=True)

# save dataframes to tab-delimited text files
gene_counts.to_csv("".join([outdir,"/",genes_output]), sep="\t", index=False)
exon_counts.to_csv("".join([outdir,"/",exons_output]), sep="\t", index=False)
intron_counts.to_csv("".join([outdir,"/",introns_output]), sep="\t", index=False)
five_prime_counts.to_csv("".join([outdir,"/",five_prime_output]), sep="\t", index=False)
intergenic_counts.to_csv("".join([outdir,"/",intergenic_output]), sep="\t", index=False)
binned_counts.to_csv("".join([outdir,"/",binned_output]), sep="\t", index=False)
