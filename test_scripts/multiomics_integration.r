# Multi-omics Integration Script

library(GenomicRanges)
library(ggplot2)
library(ComplexHeatmap)

# Load data
hic_diff <- read.csv("uninfected_vs_6hpi_diffHic_results.csv")
rna_diff <- read.csv("rna_seq_differential_expression.csv")
cutandtag <- read.csv("cutandtag_peaks.csv")
medip <- read.csv("medip_seq_methylation.csv")

# Create GRanges objects
hic_gr <- GRanges(seqnames=hic_diff$chr1, ranges=IRanges(start=hic_diff$start1, end=hic_diff$end1))
rna_gr <- GRanges(seqnames=rna_diff$chr, ranges=IRanges(start=rna_diff$start, end=rna_diff$end))
cutandtag_gr <- GRanges(seqnames=cutandtag$chr, ranges=IRanges(start=cutandtag$start, end=cutandtag$end))
medip_gr <- GRanges(seqnames=medip$chr, ranges=IRanges(start=medip$start, end=medip$end))

# Find overlaps
hic_rna_overlaps <- findOverlaps(hic_gr, rna_gr)
hic_cutandtag_overlaps <- findOverlaps(hic_gr, cutandtag_gr)
hic_medip_overlaps <- findOverlaps(hic_gr, medip_gr)

# Combine data
integrated_data <- data.frame(
  hic_diff,
  rna_expression = NA,
  cutandtag_signal = NA,
  methylation = NA
)

integrated_data$rna_expression[queryHits(hic_rna_overlaps)] <- rna_diff$log2FoldChange[subjectHits(hic_rna_overlaps)]
integrated_data$cutandtag_signal[queryHits(hic_cutandtag_overlaps)] <- cutandtag$signal[subjectHits(hic_cutandtag_overlaps)]
integrated_data$methylation[queryHits(hic_medip_overlaps)] <- medip$methylation_level[subjectHits(hic_medip_overlaps)]

# Correlation analysis
cor_matrix <- cor(integrated_data[, c("logFC", "rna_expression", "cutandtag_signal", "methylation")], use="complete.obs")

# Visualize correlations
pdf("multi_omics_correlations.pdf", width=10, height=8)
Heatmap(cor_matrix, name="Correlation", col=colorRamp2(c(-1, 0, 1), c("blue", "white", "red")))
dev.off()

# Identify regions with concordant changes across all data types
concordant_regions <- integrated_data[
  sign(integrated_data$logFC) == sign(integrated_data$rna_expression) &
  sign(integrated_data$logFC) == sign(integrated_data$cutandtag_signal) &
  sign(integrated_data$logFC) != sign(integrated_data$methylation),
]

write.csv(concordant_regions, "concordant_multi_omics_regions.csv", row.names=FALSE)

# Visualize multi-omics profile for top concordant regions
top_regions <- head(concordant_regions[order(abs(concordant_regions$logFC), decreasing=TRUE), ], 10)

pdf("top_concordant_regions_profile.pdf", width=12, height=8)
ggplot(top_regions, aes(x=reorder(interaction(chr1, start1), logFC))) +
  geom_bar(aes(y=logFC, fill="Hi-C"), stat="identity", position="dodge") +
  geom_point(aes(y=rna_expression, color="RNA-seq")) +
  geom_point(aes(y=cutandtag_signal, color="CUT&Tag")) +
  geom_point(aes(y=methylation, color="MeDIP-seq")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(x="Genomic Region", y="Normalized Signal", title="Multi-omics Profile of Top Concordant Regions") +
  scale_fill_manual(values=c("Hi-C"="blue")) +
  scale_color_manual(values=c("RNA-seq"="red", "CUT&Tag"="green", "MeDIP-seq"="purple"))
dev.off()

print("Multi-omics integration analysis complete!")
