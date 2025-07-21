library(HiCDCPlus)
library(dplyr)
library(DESeq2)
library(pheatmap)
library(stringr)
library(RColorBrewer)

outdir <- "F:/toxo_project/HiC/Hsapien_output/output_files"

#Generate features
#construct_features(output_path = paste0("F:/organism_genome/Hsapien_autosomes/GRCh38_MboI_50kb"), 
#                   gen = "Hsapiens", 
#                   gen_ver = "hg38", 
#                   sig = "GATC", 
#                   bin_type = "Bins-uniform", 
#                   binsize = 50000, 
#                   chrs = c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", 
#                            "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", 
#                            "chr17", "chr18", "chr19", "chr20", "chr21", "chr22"))

matrix_paths <- c(paste0(outdir, "/HFF_A_HiC/hicpro_files/50000/HFF_A_HiC_50000_cpm.matrix"), 
                  paste0(outdir, "/HFF_B_HiC/hicpro_files/50000/HFF_B_HiC_50000_cpm.matrix"), 
                  #paste0(outdir, "/HFF_C_HiC/hicpro_files/50000/HFF_C_HiC_50000.matrix"), 
                  #paste0(outdir, "/ME49RFP_6hpi_A_HiC/hicpro_files/50000/ME49RFP_6hpi_A_HiC_50000.matrix"), 
                  paste0(outdir, "/ME49RFP_6hpi_B_HiC/hicpro_files/50000/ME49RFP_6hpi_B_HiC_50000_cpm.matrix"), 
                  paste0(outdir, "/ME49RFP_6hpi_C_HiC/hicpro_files/50000/ME49RFP_6hpi_C_HiC_50000_cpm.matrix"), 
                  #paste0(outdir, "/ME49RFP_24hpi_A_HiC/hicpro_files/50000/ME49RFP_24hpi_A_HiC_50000.matrix"), 
                  paste0(outdir, "/ME49RFP_24hpi_B_HiC/hicpro_files/50000/ME49RFP_24hpi_B_HiC_50000_cpm.matrix"), 
                  paste0(outdir, "/ME49RFP_24hpi_C_HiC/hicpro_files/50000/ME49RFP_24hpi_C_HiC_50000_cpm.matrix"))
indexfile <- data.frame()

for(mp in matrix_paths){
    print(mp)
    output_path <- gsub('.matrix', '.txt.gz', mp)
    absfile <- gsub('_cpm.matrix', '_abs.bed', mp)
    gi_list <- generate_bintolen_gi_list(bintolen_path = "F:/organism_genome/Hsapien_autosomes/GRCh38_MboI_50kb_bintolen.txt.gz", 
                                         gen = "Hsapiens", 
                                         gen_ver = "hg38")
    gi_list <- add_hicpro_matrix_counts(gi_list, 
                                        absfile_path = absfile, 
                                        matrixfile_path = mp, 
                                        add_inter = FALSE)
    gi_list <- expand_1D_features(gi_list)
    set.seed(1010)
    gi_list <- HiCDCPlus(gi_list, ssize = 0.1)
    for (i in seq(length(gi_list))){
        indexfile <- unique(rbind(indexfile, as.data.frame(gi_list[[i]][gi_list[[i]]$qvalue<=0.05])[c('seqnames1', 'start1', 'start2')]))
    }
    gi_list_write(gi_list, fname = output_path)
}

colnames(indexfile) <- c('chr', 'startI', 'startJ')
data.table::fwrite(indexfile, 
                   paste0(outdir, '/analysis_indices.txt.gz'), 
                   sep = '\t', row.names = FALSE, quote = FALSE)

hicdcdiff(input_paths = list(HFF = c(paste0(outdir, '/HFF_A_HiC/hicpro_files/50000/HFF_A_HiC_50000.txt.gz'), 
                                     paste0(outdir, '/HFF_B_HiC/hicpro_files/50000/HFF_B_HiC_50000.txt.gz')), 
                                     #paste0(outdir, '/HFF_C_HiC/hicpro_files/50000/HFF_C_HiC_50000.txt.gz')), 
                             ME49RFP_6hpi = c(#paste0(outdir, '/ME49RFP_6hpi_A_HiC/hicpro_files/50000/ME49RFP_6hpi_A_HiC_50000.txt.gz'), 
                                              paste0(outdir, '/ME49RFP_6hpi_B_HiC/hicpro_files/50000/ME49RFP_6hpi_B_HiC_50000.txt.gz'), 
                                              paste0(outdir, '/ME49RFP_6hpi_C_HiC/hicpro_files/50000/ME49RFP_6hpi_C_HiC_50000.txt.gz')), 
                             ME49RFP_24hpi = c(#paste0(outdir, '/ME49RFP_6hpi_A_HiC/hicpro_files/50000/ME49RFP_6hpi_A_HiC_50000.txt.gz'), 
                                               paste0(outdir, '/ME49RFP_6hpi_B_HiC/hicpro_files/50000/ME49RFP_6hpi_B_HiC_50000.txt.gz'), 
                                               paste0(outdir, '/ME49RFP_6hpi_C_HiC/hicpro_files/50000/ME49RFP_6hpi_C_HiC_50000.txt.gz'))), 
          filter_file = paste0(outdir, '/analysis_indices.txt.gz'),
          output_path = paste0(outdir, '/diff_analysis/'),
          fitType = 'local',
          binsize = 50000,
          diagnostics = TRUE)





for(mp in matrix_paths){
    tads <- gi_list_topdom(gi_list, window.size = 5, file_out=true, fpath=paste0(outdir, '/diff_analysis/'))

    extract_hic_eigenvectors(
        hicfile = system.file("extdata", "eigenvector_example.hic", package = "HiCDCPlus"), 
        mode = "KR", 
        binsize = 50000, 
        gen = "Hsapiens", 
        gen_ver = "hg38", 
    )
}





source("F:/scripts/hic_scripts/HICDCPlus_Pfalciparum.R")

#Generate features
#construct_features_v2(output_path = paste0("F:/organism_genome/Pfalciparum/Pfalciparum_10kb_GATC"), 
#                      gen = "Pfalciparum", 
#                      gen_ver = "v24", 
#                      sig = c("GATC","GANTC"), 
#                      bin_type = "Bins-uniform",
#                      binsize = 10000)

#Add .hic counts
matrix_paths<-c("F:/Deitsch_collab/HiC/matrices/merged/A3/A3_10000.matrix",
                "F:/Deitsch_collab/HiC/matrices/merged/B10/B10_10000.matrix",
                "F:/Deitsch_collab/HiC/matrices/merged/D2/D2_10000.matrix",
                "F:/Deitsch_collab/HiC/matrices/merged/G91/G91_10000.matrix")
indexfile <- data.frame()

for(mp in matrix_paths){
    output_path <- paste0(gsub("^(.*[\\/])", "", gsub('.matrix', '.txt.gz', mp)))
    gi_list <- generate_bintolen_gi_list_v2(bintolen_path="F:/organism_genome/Pfalciparum3D7/Pfalciparum_10kb_GATC_bintolen.txt", gen = "Pfalciparum", gen_ver = "v24")
    gi_list <- add_hicpro_matrix_counts(gi_list, absfile_path = "F:/Deitsch_collab/HiC/matrices/merged/A3/A3_10000_abs.bed", matrixfile_path = mp, add_inter = TRUE)
    gi_list <- expand_1D_features(gi_list)
    set.seed(1010)
    gi_list <- HiCDCPlus(gi_list, ssize = 0.1)
    for (i in seq(length(gi_list))){
        indexfile <- unique(rbind(indexfile, as.data.frame(gi_list[[i]][gi_list[[i]]$qvalue<=0.05])[c('seqnames1','start1','start2')]))
    }
    gi_list_write(gi_list, fname = output_path)
}

for(hicfile_path in hicfile_paths){
    output_path <- paste0(gsub("^(.*[\\/])", "", gsub('.hic', '.txt.gz', hicfile_path)))
    
    #generate gi_list instance
    gi_list <- generate_bintolen_gi_list_v2(bintolen_path="F:/organism_genome/Pfalciparum3D7/Pfalciparum_10kb_GATC_bintolen.txt",gen="Pfalciparum",gen_ver="v24")
    gi_list <- add_hic_counts(gi_list, hic_path = hicfile_path, add_inter = TRUE)
    
    #expand features for modeling
    gi_list <- expand_1D_features(gi_list)
    set.seed(1010)
    
    #HiC-DC downsamples rows for modeling
    gi_list<-HiCDCPlus(gi_list, ssize = 0.1)
    for (i in seq(length(gi_list))){
        indexfile <- unique(rbind(indexfile, as.data.frame(gi_list[[i]][gi_list[[i]]$qvalue<=0.05])[c('seqnames1','start1','start2')]))
    }
    
    #write results to a text file
    gi_list_write(gi_list, fname = output_path)
}

#Save index file
colnames(indexfile) <- c('chr','startI','startJ')
data.table::fwrite(indexfile, "analysis_indices.txt.gz", sep='\t', row.names = FALSE, quote = FALSE)

#Differential analysis
hicdcdiff(input_paths = list(NF54 = c("nf.allValidPairs.txt.gz", "NF54A.allValidPairs.txt.gz", "NF54B.allValidPairs.txt.gz"),
                             VAR2CSA = c("var2csa.allValidPairs.txt.gz", "CSAA.allValidPairs.txt.gz", "CSAB.allValidPairs.txt.gz")),
          filter_file = "analysis_indices.txt.gz",
          output_path = "diff_analysis/",
          fitType = 'mean',
          binsize = 10000,
          diagnostics = TRUE)

##################################################################################################################################
################################################ DIFFERENTIAL INTERACTION HEATMAP ################################################
##################################################################################################################################

sizes <- read.csv('F:/organism_genome/Pfalciparum3D7/PlasmoDB-50_Pfalciparum3D7_Genome_noAPIorMIT.chrom.sizes', header = FALSE, sep = "\t")[1:14,]
colnames(sizes) <- c("chr", "length")
starts <- data.frame(matrix(NA, nrow = nrow(sizes), ncol = ncol(sizes)))
colnames(starts) <- c("chr", "start_bin")
for (i in 1:nrow(sizes)) {
  starts[i, 1] <- sizes[i, 1]
  starts[1, 2] <- 0
  if (i > nrow(sizes) - 1) {
    break
  } else {
    starts[i+1, 2] <- ceiling(sizes[i, 2] / 10000) + starts[i, 2]
  }
}
sizes$length <- ceiling(sizes$length / 10000)

var <- read.csv('F:/organism_genome/Pfalciparum3D7/var_genes.gff', header = FALSE, sep = "\t")[,c(1,4,5,9)]
colnames(var) <- c("chr", "start", "end", "gene")
for(i in 1:nrow(var)) {
  var[i, 2] <- ceiling(var[i, 2] / 10000)
  var[i, 3] <- ceiling(var[i, 3] / 10000)
}
var$gene <- str_extract_all(var$gene, "(?<==).*(?=;)")
colnames(var)[colnames(var) == "start"] <- "startI"

quantile_breaks <- function(xs, n = 20) {
  breaks <- quantile(xs, probs = seq(0, 1, length.out = n))
  breaks[!duplicated(breaks)]
}

pdf("heatmaps/diff_intra_heatmaps.pdf")

for (i in c("01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12", "13", "14")) {
  df <- read.csv(paste0("diff_analysis/diff_resD2overA3_Pf3D7_",i,"_v3.txt"), header = TRUE, sep = "\t")
  assign(paste0("VAR2CSAoverNF54_chr",i), df)
  df$startI <- df$startI / 10000 + 1
  df$startJ <- df$startJ / 10000 + 1
  df$startIplus1 <- df$startI + 1
  df$startIplus2 <- df$startI + 2
  df$startJplus2 <- df$startJ + 2
  df$startJplus1 <- df$startJ + 1
  df <- df[!(df$startI == df$startJ),]
  df <- df[!(df$startIplus1 == df$startJ),]
  df <- df[!(df$startIplus2 == df$startJ),]
  df <- df[!(df$startI == df$startJplus1),]
  df <- df[!(df$startI == df$startJplus2),]
  df <- df[!(df$padj > 0.05),]
  
  min_max <- max(abs(df$log2FoldChange))
  current_chrom <- unique(df$chr)
  list2env(split(df, df$chr), envir = .GlobalEnv)
  df %>% drop_na(chr)

  chr_mat <- matrix(0, sizes$length[which(sizes$chr == current_chrom)], sizes$length[which(sizes$chr == current_chrom)])
  chr_mat[as.matrix(df[c("startI", "startJ")])] <- df[["log2FoldChange"]]
  chr_lmat <- t(chr_mat)[,-ncol(chr_mat)]
  chr_mat[lower.tri(chr_mat)] <- chr_lmat[lower.tri(chr_lmat, diag = FALSE)]
  chr_mat <- apply(t(chr_mat),2,rev)
  
  rownames(chr_mat) <- 1:nrow(chr_mat)
  colnames(chr_mat) <- 1:ncol(chr_mat)
  
  j <- 1
  for (g in var[which(var$chr == current_chrom),]$startI) {
    rownames(chr_mat)[g] <- var[which(var$chr == current_chrom),]$gene[j]
    colnames(chr_mat)[g] <- var[which(var$chr == current_chrom),]$gene[j]
    j = j + 1
  }
  
  rownames(chr_mat) <- rev(rownames(chr_mat))
  
  mat_breaks <- quantile_breaks(df$log2FoldChange, n = 100)
  
  if (length(mat_breaks[which(mat_breaks < 0)]) > length(mat_breaks[which(mat_breaks > 0)])) {
    mat_breaks <- mat_breaks[which(mat_breaks < 0)]
    mat_breaks <- c(-max(abs(df$log2FoldChange)), mat_breaks, 0, -rev(mat_breaks), max(abs(df$log2FoldChange)))
  } else {
    mat_breaks <- mat_breaks[which(mat_breaks > 0)]
    mat_breaks <- c(-max(abs(df$log2FoldChange)), -rev(mat_breaks), 0, mat_breaks, max(abs(df$log2FoldChange)))
  }
  
  if (max(abs(df$log2FoldChange)) == max(abs(mat_breaks))) {
    mat_breaks <- head(mat_breaks, -1)
    mat_breaks <- tail(mat_breaks, -1)
  }
  mat_color <- colorRampPalette(c("blue", "white", "red"))(length(mat_breaks))
  
  annoRow <- data.frame(varGenes = factor(rownames(chr_mat) %in% var$gene))
  rownames(annoRow) <- rownames(chr_mat)
  ann_colors = list(varGenes = c("FALSE"="white", "TRUE"="green"))
  
  pheatmap(chr_mat, 
           main = paste0(gsub("_df", "", i), " ", k), 
           show_colnames = FALSE, 
           show_rownames = FALSE, 
           cluster_cols = FALSE, 
           cluster_rows = FALSE, 
           annotation_row = annoRow, 
           annotation_col = annoRow, 
           annotation_colors = ann_colors, 
           color = mat_color, 
           breaks = mat_breaks, 
           cellwidth = 360 / ncol(chr_mat), 
           cellheight = 360 / ncol(chr_mat), 
           border_color = FALSE)
}
dev.off()