library(ChIPQC)
library(GenomicFeatures)
library(GenomeInfoDb)
library(BiocParallel)
library(DiffBind)
library(tidyverse)
library(parallel)
library(pheatmap)
library(RColorBrewer)
register(SerialParam())


##################################################################################################
######################################## Reading peaksets ########################################
##################################################################################################

samples <- read.csv('F:/toxo_project/MeDIPseq/output/peak_calling/medip_metadata.txt', sep='\t')
medip_df <- dba(sampleSheet=samples)
medip_df <- dba.count(medip_df, bUseSummarizeOverlaps=TRUE)

##################################################################################################
############################ Establishing a model design and contrast ############################
##################################################################################################

medip_df <- dba.contrast(medip_df, reorderMeta=list(Condition="HFF"), minMembers=2)

##################################################################################################
############################## Performing the differential analysis ##############################
##################################################################################################

medip_df <- dba.analyze(medip_df)

##################################################################################################
########################### Retrieving the differentially bound sites ############################
##################################################################################################

medip_df.report <- dba.report(medip_df)
out <- as.data.frame(medip_df.report)
write.table(out, file="F:/toxo_project/MeDIPseq/output/peak_calling/diff_bind_peaks.txt", sep="\t", quote=F, row.names=F)

##################################################################################################
############################################ Plotting ############################################
##################################################################################################

pdf("F:/toxo_project/MeDIPseq/output/peak_calling/diffbind_figures.pdf")

# Correlation Heatmap
corrvals <- dba.plotHeatmap(medip_df, colScheme="Reds")

# PCA plot
dba.plotPCA(medip_df, contrast=1)

# Venn diagram
dba.plotVenn(medip_df, contrast=1)

# MA plot
dba.plotMA(medip_df, dotSize=1, bSmooth=FALSE)
dba.plotMA(medip_df, dotSize=1, bSmooth=FALSE, bXY=TRUE)

# Boxplot
dba.plotBox(medip_df)

# Volcano plot
dba.plotVolcano(medip_df)

# Heatmaps
hmap <- colorRampPalette(c("blue", "white", "red"))(n = 15)
readscores <- dba.plotHeatmap(medip_df, contrast=1, correlations=FALSE, scale="row", colScheme=hmap)

dev.off()
