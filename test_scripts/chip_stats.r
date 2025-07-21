options(repos = c(CRAN = "https://cloud.r-project.org/"))
library(rtracklayer)
library(ggplot2)
library(data.table)
library(DESeq2)
#install.packages('gridExtra')
library(gridExtra)
library(dplyr)

# Function to read broadPeak files
read_peaks <- function(file) {
  peaks <- read.table(file, header=FALSE, col.names=c("chrom", "start", "end", "name", "score", "strand", "signalValue", "pValue", "qValue"))
  GRanges(seqnames=peaks$chrom, ranges=IRanges(peaks$start, peaks$end), strand="*", mcols=peaks[,4:9])
}

# Read peak files
ko_peaks <- read_peaks("D2_H3K9me3_peaks.broadPeak")
wt_peaks <- read_peaks("A3_H3K9me3_peaks.broadPeak")

# Calculate widths
ko_widths <- width(ko_peaks)
wt_widths <- width(wt_peaks)

# Perform statistical test on peak widths
wilcox_test <- wilcox.test(ko_widths, wt_widths)

# Create peak width plot
peak_data <- data.frame(
  width = c(ko_widths, wt_widths),
  condition = rep(c("D2", "A3"), c(length(ko_widths), length(wt_widths)))
)

p1 <- ggplot(peak_data, aes(x = condition, y = width)) +
  geom_boxplot() +
  theme_minimal() +
  labs(title = "H3K9me3 Peak Widths", y = "Peak Width (bp)")

# Read coverage data
ko_cov <- fread("D2_H3K9_coverage.bedgraph")
wt_cov <- fread("A3_H3K9_coverage.bedgraph")

# Calculate cumulative coverage
calc_cumulative <- function(cov) {
  cov[, cum_cov := cumsum(V4) / sum(V4)]
  cov[, cum_dist := cumsum(V3 - V2) / sum(V3 - V2)]
  return(cov[, .(cum_dist, cum_cov)])
}

#ko_cum <- calc_cumulative(ko_cov)
#wt_cum <- calc_cumulative(wt_cov)

# Create cumulative distribution plot
#p2 <- ggplot() +
#  geom_line(data = ko_cum, aes(x = cum_dist, y = cum_cov, color = "KO")) +
#  geom_line(data = wt_cum, aes(x = cum_dist, y = cum_cov, color = "WT")) +
#  theme_minimal() +
#  labs(x = "Cumulative genomic distance", y = "Cumulative H3K9me3 coverage",
#       title = "Cumulative H3K9me3 coverage")

# Kolmogorov-Smirnov test
#ks_test <- ks.test(ko_cum$cum_cov, wt_cum$cum_cov)

# Differential enrichment analysis
# Assuming you have count data for H3K9me3 peaks
count_data <- read.csv("H3K9me3_counts2.txt", sep = "\t", row.names = 1)
coldata <- data.frame(condition = factor(c("D2", "D2", "A3", "A3")))

dds <- DESeqDataSetFromMatrix(countData = count_data,
                              colData = coldata,
                              design = ~ condition)

dds <- DESeq(dds)
res <- results(dds)

# Create differential enrichment plot
p3 <- ggplot(as.data.frame(res), aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point() +
  theme_minimal() +
  labs(title = "Differential H3K9me3 enrichment",
       x = "Log2 Fold Change (KO/WT)",
       y = "-log10(adjusted p-value)")

# Combine plots
pdf("H3K9me3_analysis_plots.pdf", width = 12, height = 12)
grid.arrange(p1, p2, p3, ncol = 2)
dev.off()

# Calculate statistics
median_ko <- median(ko_widths)
median_wt <- median(wt_widths)
wilcox_pvalue <- wilcox_test$p.value

total_genes <- nrow(res)
increased <- sum(res$log2FoldChange > 0 & res$padj < 0.05, na.rm = TRUE)
decreased <- sum(res$log2FoldChange < 0 & res$padj < 0.05, na.rm = TRUE)
percent_increased <- (increased / total_genes) * 100
percent_decreased <- (decreased / total_genes) * 100

ks_pvalue <- ks_test$p.value

# Write statistics to file
stats_output <- c(
  paste("Median KO peak width:", median_ko),
  paste("Median WT peak width:", median_wt),
  paste("Wilcoxon test p-value:", wilcox_pvalue),
  paste("Percent var/rifin genes with increased H3K9me3:", percent_increased),
  paste("Percent var/rifin genes with decreased H3K9me3:", percent_decreased),
  paste("K-S test p-value:", ks_pvalue),
  paste("Total peaks with significant increase:", increased),
  paste("Total peaks with significant decrease:", decreased)
)

writeLines(stats_output, "H3K9me3_analysis_stats.txt")

print("Analysis complete. Check 'H3K9me3_analysis_plots.pdf' for plots and 'H3K9me3_analysis_stats.txt' for statistics.")
