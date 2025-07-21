library(hicrep)
library(pheatmap)
library(data.table)
library(tidyr)
library(dplyr)
library(ComplexHeatmap)

bed1 <- read.delim("F:/toxo_project/HiC/Hsapien_output/output_files/HFF_A_HiC/hicpro_files/1000000/HFF_A_HiC_1000000.matrix", header=FALSE)
bed1 <- as.matrix(bed1)
mat1 <- bed2mat(bed1, resol="NONE")

bed2 <- read.delim("F:/toxo_project/HiC/Hsapien_output/output_files/HFF_B_HiC/hicpro_files/1000000/HFF_B_HiC_1000000.matrix", header=FALSE)
bed2 <- as.matrix(bed2)
mat2 <- bed2mat(bed2, resol="NONE")

bed3 <- read.delim("F:/toxo_project/HiC/Hsapien_output/output_files/HFF_C_HiC/hicpro_files/1000000/HFF_C_HiC_1000000.matrix", header=FALSE)
bed3 <- as.matrix(bed3)
mat3 <- bed2mat(bed3, resol="NONE")

bed4 <- read.delim("F:/toxo_project/HiC/Hsapien_output/output_files/ME49RFP_6hpi_A_HiC/hicpro_files/1000000/ME49RFP_6hpi_A_HiC_1000000.matrix", header=FALSE)
bed4 <- as.matrix(bed4)
mat4 <- bed2mat(bed4, resol="NONE")

bed5 <- read.delim("F:/toxo_project/HiC/Hsapien_output/output_files/ME49RFP_6hpi_B_HiC/hicpro_files/1000000/ME49RFP_6hpi_B_HiC_1000000.matrix", header=FALSE)
bed5 <- as.matrix(bed5)
mat5 <- bed2mat(bed5, resol="NONE")

bed6 <- read.delim("F:/toxo_project/HiC/Hsapien_output/output_files/ME49RFP_6hpi_C_HiC/hicpro_files/1000000/ME49RFP_6hpi_C_HiC_1000000.matrix", header=FALSE)
bed6 <- as.matrix(bed6)
mat6 <- bed2mat(bed6, resol="NONE")

bed7 <- read.delim("F:/toxo_project/HiC/Hsapien_output/output_files/ME49RFP_24hpi_A_HiC/hicpro_files/1000000/ME49RFP_24hpi_A_HiC_1000000.matrix", header=FALSE)
bed7 <- as.matrix(bed7)
mat7 <- bed2mat(bed7, resol="NONE")

bed8 <- read.delim("F:/toxo_project/HiC/Hsapien_output/output_files/ME49RFP_24hpi_B_HiC/hicpro_files/1000000/ME49RFP_24hpi_B_HiC_1000000.matrix", header=FALSE)
bed8 <- as.matrix(bed8)
mat8 <- bed2mat(bed8, resol="NONE")

bed9 <- read.delim("F:/toxo_project/HiC/Hsapien_output/output_files/ME49RFP_24hpi_C_HiC/hicpro_files/1000000/ME49RFP_24hpi_C_HiC_1000000.matrix", header=FALSE)
bed9 <- as.matrix(bed9)
mat9 <- bed2mat(bed9, resol="NONE")

h_value <- htrain(mat1, mat2, resol = 1000000, ubr = 50000000, range = 0:10)

HFF_A_HFF_A = get.scc(mat1, mat1, resol = 1000000, h = 1, lbr = 0, ubr = 50000000)
HFF_A_HFF_B = get.scc(mat1, mat2, resol = 1000000, h = 1, lbr = 0, ubr = 50000000)
HFF_A_HFF_C = get.scc(mat1, mat3, resol = 1000000, h = 1, lbr = 0, ubr = 50000000)
HFF_A_ME49RFP_6hpi_A = get.scc(mat1, mat4, resol = 1000000, h = 1, lbr = 0, ubr = 50000000)
HFF_A_ME49RFP_6hpi_B = get.scc(mat1, mat5, resol = 1000000, h = 1, lbr = 0, ubr = 50000000)
HFF_A_ME49RFP_6hpi_C = get.scc(mat1, mat6, resol = 1000000, h = 1, lbr = 0, ubr = 50000000)
HFF_A_ME49RFP_24hpi_A = get.scc(mat1, mat7, resol = 1000000, h = 1, lbr = 0, ubr = 50000000)
HFF_A_ME49RFP_24hpi_B = get.scc(mat1, mat8, resol = 1000000, h = 1, lbr = 0, ubr = 50000000)
HFF_A_ME49RFP_24hpi_C = get.scc(mat1, mat9, resol = 1000000, h = 1, lbr = 0, ubr = 50000000)

HFF_B_HFF_A = get.scc(mat2, mat1, resol = 1000000, h = 1, lbr = 0, ubr = 50000000)
HFF_B_HFF_B = get.scc(mat2, mat2, resol = 1000000, h = 1, lbr = 0, ubr = 50000000)
HFF_B_HFF_C = get.scc(mat2, mat3, resol = 1000000, h = 1, lbr = 0, ubr = 50000000)
HFF_B_ME49RFP_6hpi_A = get.scc(mat2, mat4, resol = 1000000, h = 1, lbr = 0, ubr = 50000000)
HFF_B_ME49RFP_6hpi_B = get.scc(mat2, mat5, resol = 1000000, h = 1, lbr = 0, ubr = 50000000)
HFF_B_ME49RFP_6hpi_C = get.scc(mat2, mat6, resol = 1000000, h = 1, lbr = 0, ubr = 50000000)
HFF_B_ME49RFP_24hpi_A = get.scc(mat2, mat7, resol = 1000000, h = 1, lbr = 0, ubr = 50000000)
HFF_B_ME49RFP_24hpi_B = get.scc(mat2, mat8, resol = 1000000, h = 1, lbr = 0, ubr = 50000000)
HFF_B_ME49RFP_24hpi_C = get.scc(mat2, mat9, resol = 1000000, h = 1, lbr = 0, ubr = 50000000)

HFF_C_HFF_A = get.scc(mat3, mat1, resol = 1000000, h = 1, lbr = 0, ubr = 50000000)
HFF_C_HFF_B = get.scc(mat3, mat2, resol = 1000000, h = 1, lbr = 0, ubr = 50000000)
HFF_C_HFF_C = get.scc(mat3, mat3, resol = 1000000, h = 1, lbr = 0, ubr = 50000000)
HFF_C_ME49RFP_6hpi_A = get.scc(mat3, mat4, resol = 1000000, h = 1, lbr = 0, ubr = 50000000)
HFF_C_ME49RFP_6hpi_B = get.scc(mat3, mat5, resol = 1000000, h = 1, lbr = 0, ubr = 50000000)
HFF_C_ME49RFP_6hpi_C = get.scc(mat3, mat6, resol = 1000000, h = 1, lbr = 0, ubr = 50000000)
HFF_C_ME49RFP_24hpi_A = get.scc(mat3, mat7, resol = 1000000, h = 1, lbr = 0, ubr = 50000000)
HFF_C_ME49RFP_24hpi_B = get.scc(mat3, mat8, resol = 1000000, h = 1, lbr = 0, ubr = 50000000)
HFF_C_ME49RFP_24hpi_C = get.scc(mat3, mat9, resol = 1000000, h = 1, lbr = 0, ubr = 50000000)

ME49RFP_6hpi_A_HFF_A = get.scc(mat4, mat1, resol = 1000000, h = 1, lbr = 0, ubr = 50000000)
ME49RFP_6hpi_A_HFF_B = get.scc(mat4, mat2, resol = 1000000, h = 1, lbr = 0, ubr = 50000000)
ME49RFP_6hpi_A_HFF_C = get.scc(mat4, mat3, resol = 1000000, h = 1, lbr = 0, ubr = 50000000)
ME49RFP_6hpi_A_ME49RFP_6hpi_A = get.scc(mat4, mat4, resol = 1000000, h = 1, lbr = 0, ubr = 50000000)
ME49RFP_6hpi_A_ME49RFP_6hpi_B = get.scc(mat4, mat5, resol = 1000000, h = 1, lbr = 0, ubr = 50000000)
ME49RFP_6hpi_A_ME49RFP_6hpi_C = get.scc(mat4, mat6, resol = 1000000, h = 1, lbr = 0, ubr = 50000000)
ME49RFP_6hpi_A_ME49RFP_24hpi_A = get.scc(mat4, mat7, resol = 1000000, h = 1, lbr = 0, ubr = 50000000)
ME49RFP_6hpi_A_ME49RFP_24hpi_B = get.scc(mat4, mat8, resol = 1000000, h = 1, lbr = 0, ubr = 50000000)
ME49RFP_6hpi_A_ME49RFP_24hpi_C = get.scc(mat4, mat9, resol = 1000000, h = 1, lbr = 0, ubr = 50000000)

ME49RFP_6hpi_A_HFF_A = get.scc(mat5, mat1, resol = 1000000, h = 1, lbr = 0, ubr = 50000000)
ME49RFP_6hpi_A_HFF_B = get.scc(mat5, mat2, resol = 1000000, h = 1, lbr = 0, ubr = 50000000)
ME49RFP_6hpi_A_HFF_C = get.scc(mat5, mat3, resol = 1000000, h = 1, lbr = 0, ubr = 50000000)
ME49RFP_6hpi_A_ME49RFP_6hpi_A = get.scc(mat5, mat4, resol = 1000000, h = 1, lbr = 0, ubr = 50000000)
ME49RFP_6hpi_A_ME49RFP_6hpi_B = get.scc(mat5, mat5, resol = 1000000, h = 1, lbr = 0, ubr = 50000000)
ME49RFP_6hpi_A_ME49RFP_6hpi_C = get.scc(mat5, mat6, resol = 1000000, h = 1, lbr = 0, ubr = 50000000)
ME49RFP_6hpi_A_ME49RFP_24hpi_A = get.scc(mat5, mat7, resol = 1000000, h = 1, lbr = 0, ubr = 50000000)
ME49RFP_6hpi_A_ME49RFP_24hpi_B = get.scc(mat5, mat8, resol = 1000000, h = 1, lbr = 0, ubr = 50000000)
ME49RFP_6hpi_A_ME49RFP_24hpi_C = get.scc(mat5, mat9, resol = 1000000, h = 1, lbr = 0, ubr = 50000000)

ME49RFP_6hpi_A_HFF_A = get.scc(mat6, mat1, resol = 1000000, h = 1, lbr = 0, ubr = 50000000)
ME49RFP_6hpi_A_HFF_B = get.scc(mat6, mat2, resol = 1000000, h = 1, lbr = 0, ubr = 50000000)
ME49RFP_6hpi_A_HFF_C = get.scc(mat6, mat3, resol = 1000000, h = 1, lbr = 0, ubr = 50000000)
ME49RFP_6hpi_A_ME49RFP_6hpi_A = get.scc(mat6, mat4, resol = 1000000, h = 1, lbr = 0, ubr = 50000000)
ME49RFP_6hpi_A_ME49RFP_6hpi_B = get.scc(mat6, mat5, resol = 1000000, h = 1, lbr = 0, ubr = 50000000)
ME49RFP_6hpi_A_ME49RFP_6hpi_C = get.scc(mat6, mat6, resol = 1000000, h = 1, lbr = 0, ubr = 50000000)
ME49RFP_6hpi_A_ME49RFP_24hpi_A = get.scc(mat6, mat7, resol = 1000000, h = 1, lbr = 0, ubr = 50000000)
ME49RFP_6hpi_A_ME49RFP_24hpi_B = get.scc(mat6, mat8, resol = 1000000, h = 1, lbr = 0, ubr = 50000000)
ME49RFP_6hpi_A_ME49RFP_24hpi_C = get.scc(mat6, mat9, resol = 1000000, h = 1, lbr = 0, ubr = 50000000)

ME49RFP_24hpi_A_HFF_A = get.scc(mat7, mat1, resol = 1000000, h = 1, lbr = 0, ubr = 50000000)
ME49RFP_24hpi_A_HFF_B = get.scc(mat7, mat2, resol = 1000000, h = 1, lbr = 0, ubr = 50000000)
ME49RFP_24hpi_A_HFF_C = get.scc(mat7, mat3, resol = 1000000, h = 1, lbr = 0, ubr = 50000000)
ME49RFP_24hpi_A_ME49RFP_6hpi_A = get.scc(mat7, mat4, resol = 1000000, h = 1, lbr = 0, ubr = 50000000)
ME49RFP_24hpi_A_ME49RFP_6hpi_B = get.scc(mat7, mat5, resol = 1000000, h = 1, lbr = 0, ubr = 50000000)
ME49RFP_24hpi_A_ME49RFP_6hpi_C = get.scc(mat7, mat6, resol = 1000000, h = 1, lbr = 0, ubr = 50000000)
ME49RFP_24hpi_A_ME49RFP_24hpi_A = get.scc(mat7, mat7, resol = 1000000, h = 1, lbr = 0, ubr = 50000000)
ME49RFP_24hpi_A_ME49RFP_24hpi_B = get.scc(mat7, mat8, resol = 1000000, h = 1, lbr = 0, ubr = 50000000)
ME49RFP_24hpi_A_ME49RFP_24hpi_C = get.scc(mat7, mat9, resol = 1000000, h = 1, lbr = 0, ubr = 50000000)

ME49RFP_24hpi_A_HFF_A = get.scc(mat8, mat1, resol = 1000000, h = 1, lbr = 0, ubr = 50000000)
ME49RFP_24hpi_A_HFF_B = get.scc(mat8, mat2, resol = 1000000, h = 1, lbr = 0, ubr = 50000000)
ME49RFP_24hpi_A_HFF_C = get.scc(mat8, mat3, resol = 1000000, h = 1, lbr = 0, ubr = 50000000)
ME49RFP_24hpi_A_ME49RFP_6hpi_A = get.scc(mat8, mat4, resol = 1000000, h = 1, lbr = 0, ubr = 50000000)
ME49RFP_24hpi_A_ME49RFP_6hpi_B = get.scc(mat8, mat5, resol = 1000000, h = 1, lbr = 0, ubr = 50000000)
ME49RFP_24hpi_A_ME49RFP_6hpi_C = get.scc(mat8, mat6, resol = 1000000, h = 1, lbr = 0, ubr = 50000000)
ME49RFP_24hpi_A_ME49RFP_24hpi_A = get.scc(mat8, mat7, resol = 1000000, h = 1, lbr = 0, ubr = 50000000)
ME49RFP_24hpi_A_ME49RFP_24hpi_B = get.scc(mat8, mat8, resol = 1000000, h = 1, lbr = 0, ubr = 50000000)
ME49RFP_24hpi_A_ME49RFP_24hpi_C = get.scc(mat8, mat9, resol = 1000000, h = 1, lbr = 0, ubr = 50000000)

ME49RFP_24hpi_A_HFF_A = get.scc(mat9, mat1, resol = 1000000, h = 1, lbr = 0, ubr = 50000000)
ME49RFP_24hpi_A_HFF_B = get.scc(mat9, mat2, resol = 1000000, h = 1, lbr = 0, ubr = 50000000)
ME49RFP_24hpi_A_HFF_C = get.scc(mat9, mat3, resol = 1000000, h = 1, lbr = 0, ubr = 50000000)
ME49RFP_24hpi_A_ME49RFP_6hpi_A = get.scc(mat9, mat4, resol = 1000000, h = 1, lbr = 0, ubr = 50000000)
ME49RFP_24hpi_A_ME49RFP_6hpi_B = get.scc(mat9, mat5, resol = 1000000, h = 1, lbr = 0, ubr = 50000000)
ME49RFP_24hpi_A_ME49RFP_6hpi_C = get.scc(mat9, mat6, resol = 1000000, h = 1, lbr = 0, ubr = 50000000)
ME49RFP_24hpi_A_ME49RFP_24hpi_A = get.scc(mat9, mat7, resol = 1000000, h = 1, lbr = 0, ubr = 50000000)
ME49RFP_24hpi_A_ME49RFP_24hpi_B = get.scc(mat9, mat8, resol = 1000000, h = 1, lbr = 0, ubr = 50000000)
ME49RFP_24hpi_A_ME49RFP_24hpi_C = get.scc(mat9, mat9, resol = 1000000, h = 1, lbr = 0, ubr = 50000000)

df <- data.frame(HFF_A=c(HFF_A_HFF_A$scc[1,1], HFF_A_HFF_B$scc[1,1], HFF_A_HFF_C$scc[1,1], 
                         HFF_A_ME49RFP_6hpi_A$scc[1,1], HFF_A_ME49RFP_6hpi_B$scc[1,1], HFF_A_ME49RFP_6hpi_C$scc[1,1],
                         HFF_A_ME49RFP_24hpi_A$scc[1,1], HFF_A_ME49RFP_24hpi_B$scc[1,1], HFF_A_ME49RFP_24hpi_C$scc[1,1]),
                 HFF_B=c(HFF_B_HFF_A$scc[1,1], HFF_B_HFF_B$scc[1,1], HFF_B_HFF_C$scc[1,1], 
                         HFF_B_ME49RFP_6hpi_A$scc[1,1], HFF_B_ME49RFP_6hpi_B$scc[1,1], HFF_B_ME49RFP_6hpi_C$scc[1,1],
                         HFF_B_ME49RFP_24hpi_A$scc[1,1], HFF_B_ME49RFP_24hpi_B$scc[1,1], HFF_B_ME49RFP_24hpi_C$scc[1,1]),
                 HFF_C=c(HFF_C_HFF_A$scc[1,1], HFF_C_HFF_B$scc[1,1], HFF_C_HFF_C$scc[1,1], 
                         HFF_C_ME49RFP_6hpi_A$scc[1,1], HFF_C_ME49RFP_6hpi_B$scc[1,1], HFF_C_ME49RFP_6hpi_C$scc[1,1],
                         HFF_C_ME49RFP_24hpi_A$scc[1,1], HFF_C_ME49RFP_24hpi_B$scc[1,1], HFF_C_ME49RFP_24hpi_C$scc[1,1]),
                 ME49RFP_6hpi_A=c(ME49RFP_6hpi_A_HFF_A$scc[1,1], ME49RFP_6hpi_A_HFF_B$scc[1,1], ME49RFP_6hpi_A_HFF_C$scc[1,1], 
                                  ME49RFP_6hpi_A_ME49RFP_6hpi_A$scc[1,1], ME49RFP_6hpi_A_ME49RFP_6hpi_B$scc[1,1], ME49RFP_6hpi_A_ME49RFP_6hpi_C$scc[1,1],
                                  ME49RFP_6hpi_A_ME49RFP_24hpi_A$scc[1,1], ME49RFP_6hpi_A_ME49RFP_24hpi_B$scc[1,1], ME49RFP_6hpi_A_ME49RFP_24hpi_C$scc[1,1]),
                 ME49RFP_6hpi_B=c(ME49RFP_6hpi_B_HFF_A$scc[1,1], ME49RFP_6hpi_B_HFF_B$scc[1,1], ME49RFP_6hpi_B_HFF_C$scc[1,1], 
                                  ME49RFP_6hpi_B_ME49RFP_6hpi_A$scc[1,1], ME49RFP_6hpi_B_ME49RFP_6hpi_B$scc[1,1], ME49RFP_6hpi_B_ME49RFP_6hpi_C$scc[1,1],
                                  ME49RFP_6hpi_B_ME49RFP_24hpi_A$scc[1,1], ME49RFP_6hpi_B_ME49RFP_24hpi_B$scc[1,1], ME49RFP_6hpi_B_ME49RFP_24hpi_C$scc[1,1]),
                 ME49RFP_6hpi_C=c(ME49RFP_6hpi_C_HFF_A$scc[1,1], ME49RFP_6hpi_C_HFF_B$scc[1,1], ME49RFP_6hpi_C_HFF_C$scc[1,1], 
                                  ME49RFP_6hpi_C_ME49RFP_6hpi_A$scc[1,1], ME49RFP_6hpi_C_ME49RFP_6hpi_B$scc[1,1], ME49RFP_6hpi_C_ME49RFP_6hpi_C$scc[1,1],
                                  ME49RFP_6hpi_C_ME49RFP_24hpi_A$scc[1,1], ME49RFP_6hpi_C_ME49RFP_24hpi_B$scc[1,1], ME49RFP_6hpi_C_ME49RFP_24hpi_C$scc[1,1]),
                 ME49RFP_24hpi_A=c(ME49RFP_24hpi_A_HFF_A$scc[1,1], ME49RFP_24hpi_A_HFF_B$scc[1,1], ME49RFP_24hpi_A_HFF_C$scc[1,1], 
                                   ME49RFP_24hpi_A_ME49RFP_6hpi_A$scc[1,1], ME49RFP_24hpi_A_ME49RFP_6hpi_B$scc[1,1], ME49RFP_24hpi_A_ME49RFP_6hpi_C$scc[1,1],
                                   ME49RFP_24hpi_A_ME49RFP_24hpi_A$scc[1,1], ME49RFP_24hpi_A_ME49RFP_24hpi_B$scc[1,1], ME49RFP_24hpi_A_ME49RFP_24hpi_C$scc[1,1]),
                 ME49RFP_24hpi_B=c(ME49RFP_24hpi_B_HFF_A$scc[1,1], ME49RFP_24hpi_B_HFF_B$scc[1,1], ME49RFP_24hpi_B_HFF_C$scc[1,1], 
                                   ME49RFP_24hpi_B_ME49RFP_6hpi_A$scc[1,1], ME49RFP_24hpi_B_ME49RFP_6hpi_B$scc[1,1], ME49RFP_24hpi_B_ME49RFP_6hpi_C$scc[1,1],
                                   ME49RFP_24hpi_B_ME49RFP_24hpi_A$scc[1,1], ME49RFP_24hpi_B_ME49RFP_24hpi_B$scc[1,1], ME49RFP_24hpi_B_ME49RFP_24hpi_C$scc[1,1]),
                 ME49RFP_24hpi_C=c(ME49RFP_24hpi_C_HFF_A$scc[1,1], ME49RFP_24hpi_C_HFF_B$scc[1,1], ME49RFP_24hpi_C_HFF_C$scc[1,1], 
                                   ME49RFP_24hpi_C_ME49RFP_6hpi_A$scc[1,1], ME49RFP_24hpi_C_ME49RFP_6hpi_B$scc[1,1], ME49RFP_24hpi_C_ME49RFP_6hpi_C$scc[1,1],
                                   ME49RFP_24hpi_C_ME49RFP_24hpi_A$scc[1,1], ME49RFP_24hpi_C_ME49RFP_24hpi_B$scc[1,1], ME49RFP_24hpi_C_ME49RFP_24hpi_C$scc[1,1]))

row.names(df) <- c("HFF_A", "HFF_B", "HFF_C", "ME49RFP_6hpi_A", "ME49RFP_6hpi_B", "ME49RFP_6hpi_C", "ME49RFP_24hpi_A", "ME49RFP_24hpi_B", "ME49RFP_24hpi_C")
colnames(df) <- c("HFF_A", "HFF_B", "HFF_C", "ME49RFP_6hpi_A", "ME49RFP_6hpi_B", "ME49RFP_6hpi_C", "ME49RFP_24hpi_A", "ME49RFP_24hpi_B", "ME49RFP_24hpi_C")
df = df[,order(ncol(df):1)]
df <- as.matrix(df)

ComplexHeatmap::Heatmap(df, border = "black", rect_gp = gpar(col = "black", lwd = 1),
                        cluster_rows = FALSE, row_names_side = "left", row_names_gp = gpar(fontsize = 12), 
                        cluster_columns = FALSE, column_names_rot = 45, column_names_gp = gpar(fontsize = 12), 
                        name = "stratum-adjusted\ncorrelation\ncoefficient", 
                        cell_fun = function(j, i, x, y, width, height, fill) {
                            grid.text(sprintf("%.2f", df[i, j]), x, y, gp = gpar(fontsize = 10))
                        })
