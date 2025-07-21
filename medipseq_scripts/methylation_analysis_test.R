library(MEDIPS)
library(edgeR)
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg38)
library(rtracklayer)
library(dplyr)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(BSgenome)

# Set working directory and parameters
setwd("F:/toxo_project/MeDIPseq/output/")
genome <- "hg38"
window_size <- 300

samples <- data.frame(
    sample_id = c("HFF_5mC_I", "HFF_5mC_J",
                  "ME49RFP_6hpi_5mC_C", "ME49RFP_6hpi_5mC_E",
                  "ME49RFP_24hpi_5mC_C", "ME49RFP_24hpi_5mC_E"),
    condition = c(rep("uninfected", 2), rep("6hpi", 2), rep("24hpi", 2)),
    medip_file = paste0(c("HFF_5mC_I_MeDIP/HFF_5mC_I_MeDIP", "HFF_5mC_J_MeDIP/HFF_5mC_J_MeDIP",
                          "ME49RFP_6hpi_5mC_C_MeDIP/ME49RFP_6hpi_5mC_C_MeDIP", "ME49RFP_6hpi_5mC_E_MeDIP/ME49RFP_6hpi_5mC_E_MeDIP",
                          "ME49RFP_24hpi_5mC_C_MeDIP/ME49RFP_24hpi_5mC_C_MeDIP", "ME49RFP_24hpi_5mC_E_MeDIP/ME49RFP_24hpi_5mC_E_MeDIP"),
                        ".bam"),
    input_file = paste0(c("HFF_input_I_MeDIP/HFF_input_I_MeDIP", "HFF_input_J_MeDIP/HFF_input_J_MeDIP",
                          "ME49RFP_6hpi_input_C_MeDIP/ME49RFP_6hpi_input_C_MeDIP", "ME49RFP_6hpi_input_E_MeDIP/ME49RFP_6hpi_input_E_MeDIP",
                          "ME49RFP_24hpi_input_C_MeDIP/ME49RFP_24hpi_input_C_MeDIP", "ME49RFP_24hpi_input_E_MeDIP/ME49RFP_24hpi_input_E_MeDIP"),
                        ".bam")
)

print("Sample information:")
print(samples)

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

# Function to create MEDIPS sets chromosome by chromosome
create_medips_sets <- function(samples, window_size = 300, 
                               batch_size = 5, 
                               temp_dir = "temp_medips") {
    require(MEDIPS)
    require(parallel)
    
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

print("Creating MEDIPS sets...")
medips_sets <- create_medips_sets(samples, window_size)
