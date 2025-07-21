import sys
import csv
import re
import pandas as pd
import numpy as np
import os

"""
fp_cyt = 0
cyt = 0
tot = 0

# arg 1 is the gff file for the genome of interest
for i in csv.reader(open(sys.argv[1], "r"), delimiter="\t"):
	tot += int(i[4]) - int(i[3])
	chrom = ""
	bp = 0
	if i[2] == "protein_coding_gene":
		if i[6] == "+":
			# arg 2 is the fasta file for the genome of interest
			for seqID, line in enumerate(open(sys.argv[2], 'r')):
				if line.startswith(">"):
					chrom = line.split(" | ")[0].strip(">")
				elif i[0] == chrom:
					if int(i[3]) - 560 <= bp and int(i[3]) - 500 > bp:
						fp_start = (int(i[3]) - 500) % 60
						b_count = 0
						for b in line:
							b_count += 1
							if 60 - fp_start <= b_count:
								if str(b) == "C":
									fp_cyt += 1
						bp += 60
					elif int(i[3]) - 500 <= bp and int(i[3]) - 60 > bp:
						for b in line:
							if str(b) == "C":
								fp_cyt += 1
						bp += 60
					elif int(i[3]) - 60 <= bp and int(i[3]) > bp:
						fp_end = int(i[3]) % 60
						b_count = 0
						for b in line:
							b_count += 1
							if 60 - fp_end > b_count:
								if str(b) == "C":
									fp_cyt += 1
							elif 60 - fp_end <= b_count:
								if str(b) == "C":
									cyt += 1
						bp += 60
					elif int(i[3]) <= bp and int(i[4]) - 60 >= bp:
						for b in line:
							if str(b) == "C":
								cyt += 1
						bp += 60



					elif int(i[4]) - 60 < bp and int(i[4]) >= bp:
						end = int(i[4]) % 60
						b_count = 0
						for b in line:
							b_count += 1
							if b_count <= end:
								if str(b) == "C":
									cyt += 1
						bp += 60
					elif int(i[4]) < bp:
						bp += 60
						break
					else:
						bp += 60
		if i[6] == "-":
			for seqID, line in enumerate(open(sys.argv[2], 'r')):
				if line.startswith(">"):
					chrom = line.split(" | ")[0].strip(">")
				elif i[0] == chrom:
					if int(i[3]) - 60 <= bp and int(i[3]) >= bp:
						start = int(i[3]) % 60
						b_count = 0
						for b in line:
							b_count += 1
							if 60 - start >= b_count:
								if str(b) == "T":
									cyt += 1
						bp += 60
					elif int(i[3]) < bp and int(i[4]) - 60 >= bp:
						for b in line:
							if str(b) == "T":
								cyt += 1
						bp += 60
					elif int(i[4]) - 60 < bp and int(i[4]) >= bp:
						end = int(i[4]) % 60
						b_count = 0
						for b in line:
							b_count += 1
							if b_count <= end:
								if str(b) == "T":
									cyt += 1
						bp += 60
					elif int(i[4]) < bp:
						bp += 60
						break
					else:
						bp += 60

print("".join(["total length: ", str(tot)]))
print("".join(["cytosines: ", str(cyt)]))
"""

fp1 = 0
fp2 = 0
tp1 = 0
tp2 = 0
fp1_list = []
fp2_list = []
tp1_list = []
tp2_list = []

for i in csv.reader(open(sys.argv[1], "r"), delimiter="\t"):
	chrom = ""
	bp = 0
	if i[6] == "+":
		for seqID, line in enumerate(open(sys.argv[2], 'r')):
			if line.startswith(">"):
				chrom = line.split(" | ")[0].strip(">")
			elif i[0] == chrom:
				if int(i[3]) - 210 <= bp and int(i[3]) - 150 > bp:
					fp1_start = (int(i[3]) - 150) % 60
					if fp1_start == 0:
						fp1_start = 60
					b_count = 0
					for b in line:
						b_count += 1
						if fp1_start <= b_count:
							fp1_list.append(b)
							if str(b) == "C":
								fp1 += 1
					bp += 60
				elif int(i[3]) - 150 <= bp and int(i[3]) - 60 > bp:
					for b in line:
						fp1_list.append(b)
						if str(b) == "C":
							fp1 += 1
					bp += 60
				elif int(i[3]) - 60 <= bp and int(i[3]) > bp:
					fp1_end = int(i[3]) % 60
					if fp1_end == 0:
						fp1_end = 60
					b_count = 0
					for b in line:
						b_count += 1
						if fp1_end > b_count:
							fp1_list.append(b)
							if str(b) == "C":
								fp1 += 1
						elif fp1_end <= b_count:
							fp2_list.append(b)
							if str(b) == "C":
								fp2 += 1
					bp += 60
				elif int(i[3]) <= bp and int(i[3]) + 90 > bp:
					for b in line:
						fp2_list.append(b)
						if str(b) == "C":
							fp2 += 1
					bp += 60
				elif int(i[3]) + 90 <= bp and int(i[3]) + 150 > bp and int(i[3]) + 150 < int(i[4]) - 150:
					fp2_end = (int(i[3]) + 150) % 60
					if fp2_end == 0:
						fp2_end = 60
					b_count = 0
					for b in line:
						b_count += 1
						if fp2_end > b_count:
							fp2_list.append(b)
							if str(b) == "C":
								fp2 += 1
					bp += 60
				elif int(i[3]) + 90 <= bp and int(i[3]) + 150 > bp and int(i[3]) + 150 >= int(i[4]) - 150:
					midpoint = (int(i[4]) - int(i[3])) / 2 + int(i[3])
					fp2_end = midpoint % 60
					b_count = 0
					for b in line:
						b_count += 1
						if fp2_end > b_count:
							fp2_list.append(b)
							if str(b) == "C":
								fp2 += 1
						elif fp2_end < b_count:
							tp1_list.append(b)
							if str(b) == "C":
								tp1 += 1
					bp += 60
				elif int(i[4]) - 210 <= bp and int(i[4]) - 150 > bp and int(i[3]) + 150 < int(i[4]) - 150:
					tp1_start = (int(i[4]) - 150) % 60
					b_count = 0
					for b in line:
						b_count += 1
						if tp1_start < b_count:
							tp1_list.append(b)
							if str(b) == "C":
								tp1 += 1
					bp += 60
				elif int(i[4]) - 150 <= bp and int(i[4]) - 60 > bp:
					for b in line:
						tp1_list.append(b)
						if str(b) == "C":
							tp1 += 1
					bp += 60
				elif int(i[4]) - 60 <= bp and int(i[4]) > bp:
					tp1_end = int(i[4]) % 60
					b_count = 0
					for b in line:
						b_count += 1
						if 60 - tp1_end >= b_count:
							tp1_list.append(b)
							if str(b) == "C":
								tp1 += 1
						elif 60 - tp1_end < b_count:
							tp2_list.append(b)
							if str(b) == "C":
								tp2 += 1
					bp += 60
				elif int(i[4]) <= bp and int(i[4]) + 90 >= bp:
					for b in line:
						tp2_list.append(b)
						if str(b) == "C":
							tp2 += 1
					bp += 60
				elif int(i[4]) + 90 < bp and int(i[4]) + 150 >= bp:
					tp2_end = (int(i[4]) + 150) % 60
					b_count = 0
					for b in line:
						b_count += 1
						if 60 - tp2_end >= b_count:
							tp2_list.append(b)
							if str(b) == "C":
								tp2 += 1
					bp += 60
				elif int(i[4]) + 150 < bp:
					break
				else:
					bp += 60
"""
	if i[6] == "-":
		for seqID, line in enumerate(open(sys.argv[2], 'r')):
			if line.startswith(">"):
				chrom = line.split(" | ")[0].strip(">")
			elif i[0] == chrom:
				if int(i[3]) - 60 <= bp and int(i[3]) >= bp:
					start = int(i[3]) % 60
					b_count = 0
					for b in line:
						b_count += 1
						if 60 - start >= b_count:
							if str(b) == "T":
								cyt += 1
					bp += 60
				elif int(i[3]) < bp and int(i[4]) - 60 >= bp:
					for b in line:
						if str(b) == "T":
							cyt += 1
					bp += 60
				elif int(i[4]) - 60 < bp and int(i[4]) >= bp:
					end = int(i[4]) % 60
					b_count = 0
					for b in line:
						b_count += 1
						if b_count <= end:
							if str(b) == "T":
								cyt += 1
					bp += 60
				elif int(i[4]) < bp:
					bp += 60
					break
				else:
					bp += 60
"""
print("".join(["five prime list 1:", str(len([i for i in fp1_list if i != "\n"]))]))
print("".join(fp1_list))
print("".join(["five prime list 2:", str(len([i for i in fp2_list if i != "\n"]))]))
print("".join(fp2_list))
print("".join(["three prime list 1:", str(len([i for i in tp1_list if i != "\n"]))]))
print("".join(tp1_list))
print("".join(["three prime list 2:", str(len([i for i in tp2_list if i != "\n"]))]))
print("".join(tp2_list))
print()
print("".join(["five prime of intron: ", str(fp1)]))
print("".join(["five prime end of intron: ", str(fp2)]))
print("".join(["three prime end of intron: ", str(tp1)]))
print("".join(["three prime of intron: ", str(tp2)]))
