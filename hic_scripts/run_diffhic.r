# Load necessary libraries
if (!requireNamespace("diffHic", quietly = TRUE)) {
    BiocManager::install("diffHic")
}
library(diffHic)
library(edgeR)  # Required by diffHic for statistical analysis

# Define file paths to Hi-C data for each condition and replicate
file_paths <- list(
    uninfected = c("/mnt/f/toxo_project/HiC/Hsapien_output/HFF_A_HiC/juicer_files/HFF_A_HiC.hic", "/mnt/f/toxo_project/HiC/Hsapien_output/HFF_B_HiC/juicer_files/HFF_B_HiC.hic", "/mnt/f/toxo_project/HiC/Hsapien_output/HFF_C_HiC/juicer_files/HFF_C_HiC.hic"),
    ME49RFP_6hpi = c("/mnt/f/toxo_project/HiC/Hsapien_output/ME49RFP_6hpi_A_HiC/juicer_files/ME49RFP_6hpi_A_HiC.hic", "/mnt/f/toxo_project/HiC/Hsapien_output/ME49RFP_6hpi_B_HiC/juicer_files/ME49RFP_6hpi_B_HiC.hic", "/mnt/f/toxo_project/HiC/Hsapien_output/ME49RFP_6hpi_C_HiC/juicer_files/ME49RFP_6hpi_C_HiC.hic"),
    ME49RFP_24hpi = c("/mnt/f/toxo_project/HiC/Hsapien_output/ME49RFP_24hpi_A_HiC/juicer_files/ME49RFP_24hpi_A_HiC.hic", "/mnt/f/toxo_project/HiC/Hsapien_output/ME49RFP_24hpi_B_HiC/juicer_files/ME49RFP_24hpi_B_HiC.hic", "/mnt/f/toxo_project/HiC/Hsapien_output/ME49RFP_24hpi_C_HiC/juicer_files/ME49RFP_24hpi_C_HiC.hic")
)

# Set bin size (adjust as needed, e.g., 1 Mb = 1e6, 100 kb = 1e5)
bin_size <- 1e6  # 1 Mb binning

# Step 1: Import and bin the Hi-C data
# Load and bin each Hi-C replicate, then combine them
hic_data <- lapply(file_paths, function(files) {
    lapply(files, function(file) importContacts(file, binSize=bin_size))
})

# Combine all samples into a single DGEList object
y <- do.call(c, unlist(hic_data, recursive=FALSE))

# Step 2: Set up sample information and design matrix
# Define conditions and replicates
conditions <- factor(rep(c("uninfected", "6hpi", "24hpi"), each=3))
design <- model.matrix(~ conditions)

# Step 3: Normalize the Hi-C data
y <- normOffsets(y, method="loess")

# Step 4: Estimate dispersions
y <- estimateDisp(y, design)

# Step 5: Fit a generalized linear model (GLM) to the data
fit <- glmFit(y, design)

# Step 6: Perform differential interaction testing for time contrasts
# Define contrasts for each pairwise comparison
contrast_6h_vs_uninfected <- makeContrasts(conditions6hpi - conditionsuninfected, levels=design)
contrast_24h_vs_6h <- makeContrasts(conditions24hpi - conditions6hpi, levels=design)
contrast_24h_vs_uninfected <- makeContrasts(conditions24hpi - conditionsuninfected, levels=design)

# Perform likelihood ratio tests for each contrast
lrt_6h_vs_uninfected <- glmLRT(fit, contrast=contrast_6h_vs_uninfected)
lrt_24h_vs_6h <- glmLRT(fit, contrast=contrast_24h_vs_6h)
lrt_24h_vs_uninfected <- glmLRT(fit, contrast=contrast_24h_vs_uninfected)

# Step 7: Extract significant interactions
fdr_threshold <- 0.05

# Filter significant interactions for each comparison
significant_interactions_6h_vs_uninfected <- topTags(lrt_6h_vs_uninfected, n=Inf, p.value=fdr_threshold)
significant_interactions_24h_vs_6h <- topTags(lrt_24h_vs_6h, n=Inf, p.value=fdr_threshold)
significant_interactions_24h_vs_uninfected <- topTags(lrt_24h_vs_uninfected, n=Inf, p.value=fdr_threshold)

# Step 8: Save results to files
write.table(significant_interactions_6h_vs_uninfected$table, "differential_interactions_6h_vs_uninfected.txt", sep="\t", row.names=FALSE, quote=FALSE)
write.table(significant_interactions_24h_vs_6h$table, "differential_interactions_24h_vs_6h.txt", sep="\t", row.names=FALSE, quote=FALSE)
write.table(significant_interactions_24h_vs_uninfected$table, "differential_interactions_24h_vs_uninfected.txt", sep="\t", row.names=FALSE, quote=FALSE)

# Optional: Visualize results using MDS plot or heatmaps
plotMDS(y, col=c("red", "green", "blue")[as.integer(conditions)])
