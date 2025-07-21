import os, re, sys, csv
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

# Generates dataframes from input bed files and assigns output file 
def assign_args(args):

	# Iterate over input arguments
	for i in args:

		# Input is gff containing gene info
		if os.path.splitext(i)[1] == ".gff":

			# Generate pandas dataframes from input gff
			gff = pd.read_csv(i, sep='\t', header=None)

			#################################################################################################
			#gff = gff[gff[0] == "Pf3D7_01_v3"]
			#################################################################################################

		# Input is bed file containing read count at each nucleotide
		elif os.path.splitext(i)[1] == ".bed":

			# Generate pandas dataframe from input peaks file
			bed = pd.read_csv(i, sep='\t', header=None)

			#################################################################################################
			#bed = bed[bed[0] == "Pf3D7_01_v3"]
			#################################################################################################

			# Assign output file
			outfile = "".join([os.path.splitext(os.path.basename(i))[0], "_introns.pdf"])
			output = "".join([os.path.split(os.path.dirname(os.path.abspath(i)))[0], "/", outfile])

		# If any input is not a GFF or bed format, 
		# an error is raised and the script exits immediately.
		else:
			print()
			print("".join(["Input file ", i, " is not a recognized input format."]))
			print("Run script with standard GFF and bed files.")
			print()
			sys.exit(0)

	return gff, bed, output

# Filters gff to generate dataframe of gene coordinates and extracts intron count per gene
def filter_gff(gff):

	# Filter gff dataframe to only include genes or exons
	genes = gff[(gff[2] == "protein_coding_gene") | (gff[2] == "ncRNA_gene") | (gff[2] == "pseudogene")].copy()
	exons = gff[(gff[2] == "exon")].copy()

	# Extract accession numbers for genes and exons
	genes[9] = genes[8].str.extract(r'ID=(.*?);', expand=True).copy()
	exons[9] = exons[8].str.extract(r'ID=exon_(.*?).[1-9]-', expand=True).copy()

	# Drop unnecessary columns from gene dataframe
	genes.drop([1,2,5,6,7,8], axis=1, inplace=True)
	genes = genes.rename(columns={0:"chr", 3:"start", 4:"end", 9:"gene"})

	# Find number of introns in each gene
	intron_counts = exons[9].value_counts().reset_index().sort_values("index").rename(columns={"index":"gene", 9:"introns"}).reset_index(drop=True)
	intron_counts["introns"] -= 1

	return genes, intron_counts

def gene_counts(genes, bed):
	count_list = []
	for chrom in genes["chr"].unique():
		sub_bed = bed[bed[0] == chrom]
		sub_genes = genes[genes["chr"] == chrom]
		for gene in sub_genes.itertuples():
			gene_bed = sub_bed.loc[(sub_bed[0] == gene[1]) & (sub_bed[1] >= gene[2]) & (sub_bed[1] <= gene[3])].reset_index(drop=True).copy()
			count_list.append(gene_bed[2].sum() / len(gene_bed.index))
	return count_list

def intron_reads(introns, genes):
	introns["reads"] = genes
	count_sums = introns.groupby("introns")["reads"].sum().tolist()
	intron_nums = introns["introns"].value_counts().reset_index(drop=True).tolist()
	intron_counts = [a/b for a, b in zip(count_sums,intron_nums)]
	return intron_counts

def plot_introns(introns, output):
	pdf = PdfPages(output)
	barplt = plt.bar(np.array(list(range(len(introns)))), np.array(introns), color="#f8766d") #00bfc4
	pdf.savefig()
	plt.close()
	pdf.close()

def main():
	# Assign input and output args
	gff, bed, output= assign_args(sys.argv[1:])

	# Filter gff to get gene coordinates and intron counts
	genes, intron_counts = filter_gff(gff)

	# Find read count per gene
	gene_read_count = gene_counts(genes, bed)

	# Find read count based on intron count
	intron_read_count = intron_reads(intron_counts, gene_read_count)
	
	# Plot barplot with read counts based on intron count
	plot_introns(intron_read_count, output)

if __name__ == "__main__":
	main()