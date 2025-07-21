# MeDIP-seq Differential Methylation Analysis
# Author: Todd Lenz, todd.lenz9988@gmail.com
# Date Created: June 6, 2025
# Date Updated: June 7, 2025

library(MEDIPS)
library(edgeR)
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg38)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(rtracklayer)
library(dplyr)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(BSgenome)
library(parallel)
library(ChIPseeker)
library(clusterProfiler)

# Set working directory and parameters
setwd("F:/toxo_project/MeDIPseq/output/")
genome <- "hg38"
window_size <- 500

samples <- data.frame(
    sample_id = c("HFF_5mC_rep1", "HFF_5mC_rep2",
                  "ME49RFP_6hpi_5mC_rep1", "ME49RFP_6hpi_5mC_rep2",
                  "ME49RFP_24hpi_5mC_rep1", "ME49RFP_24hpi_5mC_rep2"),
    condition = c(rep("uninfected", 2), rep("6hpi", 2), rep("24hpi", 2)),
    medip_file = paste0(c("HFF_5mC_I_MeDIP/HFF_5mC_I_MeDIP", "HFF_5mC_J_MeDIP/HFF_5mC_J_MeDIP",
                          "ME49RFP_6hpi_5mC_C_MeDIP/ME49RFP_6hpi_5mC_C_MeDIP", "ME49RFP_6hpi_5mC_E_MeDIP/ME49RFP_6hpi_5mC_E_MeDIP_seed4",
                          "ME49RFP_24hpi_5mC_C_MeDIP/ME49RFP_24hpi_5mC_C_MeDIP", "ME49RFP_24hpi_5mC_E_MeDIP/ME49RFP_24hpi_5mC_E_MeDIP"),
                        ".bam"),
    input_file = paste0(c("HFF_input_I_MeDIP/HFF_input_I_MeDIP_seed4", "HFF_input_J_MeDIP/HFF_input_J_MeDIP_seed4",
                          "ME49RFP_6hpi_input_C_MeDIP/ME49RFP_6hpi_input_C_MeDIP", "ME49RFP_6hpi_input_E_MeDIP/ME49RFP_6hpi_input_E_MeDIP_seed4",
                          "ME49RFP_24hpi_input_C_MeDIP/ME49RFP_24hpi_input_C_MeDIP", "ME49RFP_24hpi_input_E_MeDIP/ME49RFP_24hpi_input_E_MeDIP_seed4"),
                        ".bam")
)

################################################################################
############################## Load MeDIPseq data ##############################
################################################################################

# Generate MEDIPS dataset for downstream analysis
get_chromosomes_from_bam <- function(bam_file) {
    require(Rsamtools)
    bam_header <- scanBamHeader(bam_file)
    chrom_names <- names(bam_header[[1]]$targets)
    
    standard_chroms <- chrom_names[grepl("^chr[0-9XY]+$", chrom_names)]
    return(standard_chroms)
}

# Function to create MEDIPS sets
create_medips_sets <- function(samples, window_size = 500, 
                               batch_size = 5, 
                               temp_dir = "temp_medips") {
    
    if (!dir.exists(temp_dir)) {
        dir.create(temp_dir, recursive = TRUE)
    }
    
    chromosomes <- get_chromosomes_from_bam(samples$medip_file[1])
    cat("Found chromosomes:", paste(chromosomes, collapse = ", "), "\n")
    
    chr_batches <- split(chromosomes, ceiling(seq_along(chromosomes) / batch_size))
    
    medip_sets <- list()
    input_sets <- list()
    
    for(i in 1:nrow(samples)) {
        sample_id <- samples$sample_id[i]
        medip_file <- samples$medip_file[i]
        input_file <- samples$input_file[i]
        
        cat("Processing sample:", sample_id, "\n")
        
        # Process chromosomes in batches
        medip_chr_data <- list()
        input_chr_data <- list()
        
        for(batch_idx in seq_along(chr_batches)) {
            chr_batch <- chr_batches[[batch_idx]]
            cat("Processing chromosome batch", batch_idx, ":", paste(chr_batch, collapse = ", "), "\n")
            
            medip_batch <- MEDIPS.createSet(
                file = medip_file,
                BSgenome = "BSgenome.Hsapiens.UCSC.hg38",
                extend = 0,
                shift = 0,
                window_size = window_size,
                uniq = 1e-3,
                chr.select = chr_batch,
                paired = TRUE
            )
            
            input_batch <- MEDIPS.createSet(
                file = input_file,
                BSgenome = "BSgenome.Hsapiens.UCSC.hg38",
                extend = 0,
                shift = 0,
                window_size = window_size,
                uniq = 1e-3,
                chr.select = chr_batch,
                paired = TRUE
            )
            
            medip_chr_data[[batch_idx]] <- medip_batch
            input_chr_data[[batch_idx]] <- input_batch
            
            gc(verbose = FALSE)
        }
        
        medip_sets[[sample_id]] <- medip_chr_data
        input_sets[[sample_id]] <- input_chr_data
        
        cat("  Completed sample:", sample_id, "\n")
    }
    
    return(list(medip = medip_sets, input = input_sets))
}

merge_medips_batches <- function(samples, medips_sets, window_size = 500) {
    combined_medip <- list()
    combined_input <- list()
    
    for (sample in samples$sample_id) {
        batch_size <- length(medips_sets$medip[[sample]])
        
        chr_lengths_medip <- list()
        chr_lengths_input <- list()
        chr_names_medip <- character()
        chr_names_input <- character()
        total_regions_medip <- 0
        total_regions_input <- 0
        all_counts_medip <- list()
        all_counts_input <- list()
        
        for (batch in 1:batch_size) {
            medip_batch <- medips_sets$medip[[sample]][[batch]]
            input_batch <- medips_sets$input[[sample]][[batch]]
            
            chr_lengths_medip <- c(chr_lengths_medip, slot(medip_batch, "chr_lengths"))
            chr_lengths_input <- c(chr_lengths_input, slot(input_batch, "chr_lengths"))
            
            chr_names_medip <- unique(c(chr_names_medip, slot(medip_batch, "chr_names")))
            chr_names_input <- unique(c(chr_names_input, slot(input_batch, "chr_names")))
            
            total_regions_medip <- total_regions_medip + slot(medip_batch, "number_regions")
            total_regions_input <- total_regions_input + slot(input_batch, "number_regions")
            
            all_counts_medip[[batch]] <- as.numeric(slot(medip_batch, "genome_count"))
            all_counts_input[[batch]] <- as.numeric(slot(input_batch, "genome_count"))
        }
        
        # Merge chr lengths and counts
        merged_chr_lengths_medip <- unlist(chr_lengths_medip)
        merged_chr_lengths_input <- unlist(chr_lengths_input)
        
        merged_counts_medip <- do.call(c, all_counts_medip)
        merged_counts_input <- do.call(c, all_counts_input)
        
        # Create new MEDIPSset objects
        merged_medip <- new("MEDIPSset",
                            sample_name = slot(medips_sets$medip[[sample]][[1]], "sample_name"),
                            path_name = slot(medips_sets$medip[[sample]][[1]], "path_name"),
                            genome_name = "BSgenome.Hsapiens.UCSC.hg38",
                            number_regions = total_regions_medip,
                            chr_names = chr_names_medip,
                            chr_lengths = merged_chr_lengths_medip,
                            window_size = window_size,
                            extend = 0,
                            shifted = 0,
                            uniq = 1e-3,
                            genome_count = merged_counts_medip)
        
        merged_input <- new("MEDIPSset",
                            sample_name = slot(medips_sets$input[[sample]][[1]], "sample_name"),
                            path_name = slot(medips_sets$input[[sample]][[1]], "path_name"),
                            genome_name = "BSgenome.Hsapiens.UCSC.hg38",
                            number_regions = total_regions_input,
                            chr_names = chr_names_input,
                            chr_lengths = merged_chr_lengths_input,
                            window_size = window_size,
                            extend = 0,
                            shifted = 0,
                            uniq = 1e-3,
                            genome_count = merged_counts_input)
        
        combined_medip[[sample]] <- merged_medip
        combined_input[[sample]] <- merged_input
    }
    
    return(list(medip = combined_medip, input = combined_input))
}

print("Creating MEDIPS sets...")
medips_sets <- create_medips_sets(samples, window_size)
merged_medips_sets <- merge_medips_batches(samples, medips_sets)
names(merged_medips_sets$medip) <- c("HFF_5mC_rep1", "HFF_5mC_rep2",
                                     "6hpi_5mC_rep1", "6hpi_5mC_rep2",
                                     "24hpi_5mC_rep1", "24hpi_5mC_rep2")
names(merged_medips_sets$input) <- c("HFF_input_rep1", "HFF_input_rep2",
                                     "6hpi_input_rep1", "6hpi_input_rep2",
                                     "24hpi_input_rep1", "24hpi_input_rep2")
sample_names <- c("HFF_5mC_rep1", "HFF_5mC_rep2",
                  "6hpi_5mC_rep1", "6hpi_5mC_rep2",
                  "24hpi_5mC_rep1", "24hpi_5mC_rep2")

for (i in seq_along(merged_medips_sets$medip)) {
    merged_medips_sets$medip[[i]]@sample_name <- sample_names[i]
}


################################################################################
######################### QC analysis of MeDIPseq data #########################
################################################################################

pearson_cor <- MEDIPS.correlation(MSets=merged_medips_sets$medip,
                                  method = "pearson", plot=FALSE)
pearson_cor[lower.tri(pearson_cor)] <- t(pearson_cor)[lower.tri(pearson_cor)]
colors <- colorRampPalette(brewer.pal(9, "Reds"))(255)
pheatmap(pearson_cor,
         display_numbers = TRUE,
         col = colors,
         number_color = "#000000",
         main = "Pearson Correlation",
         treeheight_row = 0,
         legend = FALSE)

spearman_cor <- MEDIPS.correlation(MSets=merged_medips_sets$medip,
                                   method = "spearman", plot=FALSE)
spearman_cor[lower.tri(spearman_cor)] <- t(spearman_cor)[lower.tri(spearman_cor)]
colors <- colorRampPalette(brewer.pal(9, "Reds"))(255)
pheatmap(spearman_cor,
         display_numbers = TRUE,
         col = colors,
         number_color = "#000000",
         main = "Spearman Correlation",
         treeheight_row = 0,
         legend = FALSE)


################################################################################
###################### Differential methylation analysis #######################
################################################################################

# Differential methylation analysis
run_differential_analysis <- function(medips_sets, samples,
                                      statistical_method = "edgeR",
                                      p_value = 0.1) {
    
    # Get sample information
    medip_sets <- medips_sets$medip
    input_sets <- medips_sets$input
    
    CS <- MEDIPS.couplingVector(
        pattern = "CG",
        refObj = merged_medips_sets$medip[[1]])
    
    # Find sample groups and replicate count
    conditions <- unique(samples$condition)
    pairwise_comparisons <- combn(conditions, 2, simplify = FALSE)
    
    # Perform comparisons for this batch
    diff_comp <- list()
    sig_diff <- list()
    sig_adj_diff <- list()
    
    for(comp in pairwise_comparisons) {
        group1 <- comp[1]
        group2 <- comp[2]
        cat("  Comparing", group1, "vs", group2, "\n")
        
        group1_ids <- samples$sample_id[which(samples$condition == group1)]
        group2_ids <- samples$sample_id[which(samples$condition == group2)]
        
        comp_name <- paste0(group1, "_vs_", group2)
        
        diff <- MEDIPS.meth(MSet1 = medip_sets[group1_ids],
                            MSet2 = medip_sets[group2_ids],
                            CSet = CS,
                            ISet1 = input_sets[group1_ids],
                            ISet2 = input_sets[group2_ids],
                            p.adj = "BH",
                            diff.method = "edgeR",
                            MeDIP = TRUE,
                            CNV = FALSE,
                            minRowSum = 0,
                            diffnorm = "tmm")
        diff_comp[[comp_name]] <- diff
        
        sig <- MEDIPS.selectSig(results = diff,
                                p.value = p_value,
                                adj = TRUE,
                                ratio = NULL,
                                bg.counts = NULL,
                                CNV = FALSE)
        sig_diff[[comp_name]] <- sig
        
        sig_adj <- MEDIPS.selectSig(results = diff,
                                    p.value = p_value,
                                    adj = TRUE,
                                    ratio = NULL,
                                    bg.counts = NULL,
                                    CNV = FALSE)
        sig_adj_diff[[comp_name]] <- sig_adj
    }
    return(list(diff = diff_comp, 
                sig = sig_diff,
                sig_adj = sig_adj_diff))
}

diff_results <- run_differential_analysis(merged_medips_sets, samples, 
                                          statistical_method = "edgeR", 
                                          p_value = 0.1)




diff <- MEDIPS.meth(MSet1 = c(merged_medips_sets$medip["HFF_5mC_rep1"], merged_medips_sets$medip["HFF_5mC_rep2"]),
                    MSet2 = c(merged_medips_sets$medip["ME49RFP_6hpi_5mC_rep1"], merged_medips_sets$medip["ME49RFP_6hpi_5mC_rep2"]),
                    CSet = CS,
                    ISet1 = c(merged_medips_sets$input["HFF_input_rep1"], merged_medips_sets$input["HFF_input_rep2"]),
                    ISet2 = c(merged_medips_sets$input["ME49RFP_6hpi_input_rep1"], merged_medips_sets$input["ME49RFP_6hpi_input_rep2"]),
                    p.adj = "BH",
                    diff.method = "edgeR",
                    MeDIP = TRUE,
                    CNV = FALSE,
                    minRowSum = 0,
                    diffnorm = "tmm")







valid_6hpi <- !is.na(diff_results$diff$uninfected_vs_6hpi$edgeR.p.value) & 
    !is.na(diff_results$diff$uninfected_vs_6hpi$edgeR.logFC) & 
    diff_results$diff$uninfected_vs_6hpi$edgeR.p.value < 0.05 & 
    diff_results$diff$uninfected_vs_6hpi$edgeR.logFC > 1
valid_24hpi <- !is.na(diff_results$diff$uninfected_vs_24hpi$edgeR.p.value) & 
    !is.na(diff_results$diff$uninfected_vs_24hpi$edgeR.logFC) & 
    diff_results$diff$uninfected_vs_24hpi$edgeR.p.value < 0.05 & 
    diff_results$diff$uninfected_vs_24hpi$edgeR.logFC > 1

gr_dmrs_6hpi <- GRanges(
    seqnames = diff_results$diff$uninfected_vs_6hpi$chr[valid_6hpi],
    ranges = IRanges(
        start = diff_results$diff$uninfected_vs_6hpi$start[valid_6hpi],
        end = diff_results$diff$uninfected_vs_6hpi$stop[valid_6hpi]
    )
)
gr_dmrs_24hpi <- GRanges(
    seqnames = diff_results$diff$uninfected_vs_24hpi$chr[valid_24hpi],
    ranges = IRanges(
        start = diff_results$diff$uninfected_vs_24hpi$start[valid_24hpi],
        end = diff_results$diff$uninfected_vs_24hpi$stop[valid_24hpi]
    )
)

# Annotate with ChIPseeker
peakAnno_6hpi <- annotatePeak(gr_dmrs_6hpi,
                              TxDb = TxDb.Hsapiens.UCSC.hg38.knownGene,
                              tssRegion = c(-2000, 2000),
                              annoDb = "org.Hs.eg.db")
peakAnno_24hpi <- annotatePeak(gr_dmrs_24hpi,
                               TxDb = TxDb.Hsapiens.UCSC.hg38.knownGene,
                               tssRegion = c(-2000, 2000),
                               annoDb = "org.Hs.eg.db")

plotAnnoPie(peakAnno)
plotDistToTSS(peakAnno)

gene_ids <- as.data.frame(peakAnno)$geneId
ego <- enrichGO(gene = gene_ids,
                OrgDb = org.Hs.eg.db,
                keyType = "ENTREZID",
                ont = "BP", pAdjustMethod = "BH",
                qvalueCutoff = 0.05, readable = TRUE)

dotplot(ego, showCategory = 20)

# Load DESeq2 results (must include gene symbol or Entrez ID)
rna_df_6 <- read.csv("F:/toxo_project/RNAseq/Hsapien_output/output_files/ME49RFP_6hpi_vs_HFF_DEgenes.txt")
rna_df_24 <- read.csv("F:/toxo_project/RNAseq/Hsapien_output/output_files/ME49RFP_24hpi_vs_HFF_DEgenes.txt")

# Match on gene ID
dmr_genes <- as.data.frame(peakAnno)
merged_df <- dmr_genes %>%
    left_join(rna_df, by = c("geneId" = "gene_id")) %>%
    filter(log2FoldChange < 0 & logFC > 1) # Hyper-methylated & downregulated

# Save immune-related candidates
write.csv(merged_df, "suppressed_immune_genes.csv", row.names = FALSE)
