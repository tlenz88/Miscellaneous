# Load libraries
library(tidyverse)
library(data.table)
library(ggplot2)

# --- Load data ---
data <- read.csv("F:/ready_to_be_compressed/loic_project/chip_vs_rna_read_counts.txt", sep = "\t")

data <- data %>% mutate(across(-Gene_ID, ~log2(. + 1)))

comparisons <- tribble(
    ~comparison_name,           ~chip_col,           ~rna_col,
    "Astro_D3_vs_RNA_D3",       "Astro_D3_H3K4",     "Astro_D3_RNA",
    "Astro_D3_vs_RNA_D5",       "Astro_D3_H3K4",     "Astro_D5_RNA",
    "Astro_D3_vs_RNA_D7",       "Astro_D3_H3K4",     "Astro_D7_RNA",
    "Astro_D5_vs_RNA_D5",       "Astro_D5_H3K4",     "Astro_D5_RNA",
    "Astro_D5_vs_RNA_D7",       "Astro_D5_H3K4",     "Astro_D7_RNA",
    "Astro_D7_vs_RNA_D7",       "Astro_D7_H3K4",     "Astro_D7_RNA",
    "HFF_D3_vs_RNA_D3",         "HFF_D3_H3K4",       "HFF_D3_RNA",
    "HFF_D3_vs_RNA_D5",         "HFF_D3_H3K4",       "HFF_D5_RNA",
    "HFF_D3_vs_RNA_D7",         "HFF_D3_H3K4",       "HFF_D7_RNA",
    "HFF_D5_vs_RNA_D5",         "HFF_D5_H3K4",       "HFF_D5_RNA",
    "HFF_D5_vs_RNA_D7",         "HFF_D5_H3K4",       "HFF_D7_RNA",
    "HFF_D7_vs_RNA_D7",         "HFF_D7_H3K4",       "HFF_D7_RNA"
)

compute_correlation <- function(df, chip_col, rna_col) {
    chip <- df[[chip_col]]
    rna <- df[[rna_col]]
    
    valid_idx <- complete.cases(chip, rna)
    chip <- chip[valid_idx]
    rna <- rna[valid_idx]
    gene_ids <- df$Gene_ID[valid_idx]
    
    cor_val <- cor(chip, rna, method = "spearman")
    
    # Rank similarity (1 = identical rank, 0 = opposite)
    ranks_chip <- rank(-chip)
    ranks_rna <- rank(-rna)
    rank_similarity <- 1 - abs(ranks_chip - ranks_rna) / length(ranks_chip)
    
    tibble(
        Gene_ID = gene_ids,
        H3K4me3 = chip,
        RNA = rna,
        Spearman = cor_val,
        RankSimilarity = rank_similarity
    )
}

results_list <- list()

for (i in seq_len(nrow(comparisons))) {
    comp <- comparisons[i, ]
    res <- compute_correlation(data, comp$chip_col, comp$rna_col)
    res$Comparison <- comp$comparison_name
    results_list[[comp$comparison_name]] <- res
}

all_results <- bind_rows(results_list)

fwrite(all_results, "chip_rna_correlation_results.tsv", sep = "\t")

summary_corrs <- comparisons %>%
    rowwise() %>%
    mutate(
        Spearman = cor(data[[chip_col]], data[[rna_col]], method = "spearman", use = "complete.obs")
    )

write_tsv(summary_corrs, "global_spearman_summary.tsv")
