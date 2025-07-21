library(org.Hs.eg.db)
library(clusterProfiler)
library(ReactomePA)
library(enrichplot)
library(ggplot2)
library(DESeq2)
library(dplyr)
library(cowplot)
library(stats)
library(openxlsx)
library(stats)
library(reshape2)

GSEA_analysis <- function(deg, title, ontology) {
    print("Starting GSEA analysis")
    entrezIDs <- mapIds(org.Hs.eg.db, 
                        keys = deg$geneID, 
                        column = "ENTREZID", 
                        keytype = "SYMBOL")
    deg$entrezID <- entrezIDs
    deg <- deg[, c("entrezID", "log2FoldChange")]
    deg <- na.omit(deg)
    geneList <- deg$log2FoldChange
    names(geneList) <- deg$entrezID
    geneList <- sort(geneList, decreasing = TRUE)
    print(paste("Number of genes in geneList:", length(geneList)))
    
    # Perform GSEA
    gsea <- gseGO(
        geneList = geneList,
        ont = ontology,
        OrgDb = "org.Hs.eg.db",
        keyType = "ENTREZID",
        minGSSize = 10,
        maxGSSize = 500,
        pvalueCutoff = 0.05,
        pAdjustMethod = "BH",
        verbose = TRUE,
        nPermSimple = 10000,
        eps = 0
    )
    
    print("GSEA completed")
    print(paste("Number of enriched gene sets:", nrow(gsea@result)))
    
    if (nrow(gsea@result) == 0) {
        print("No enriched gene sets found. Returning NULL.")
        return(NULL)
    }
    
    # Generate plots
    theme_opt <- theme_bw() + 
        theme(
            axis.text.y = element_text(size = 8, vjust = 0.5),
            axis.text.x = element_text(size = 8),
            axis.title.y = element_text(size = 8),
            axis.title.x = element_text(size = 8)
        )
    
    p1 <- dotplot(gsea, 
                  showCategory = min(10, nrow(gsea)), 
                  title = title)
    p1 <- p1 + theme_opt
    
    p2 <- ridgeplot(gsea, 
                    showCategory = min(10, nrow(gsea)), 
                    fill = "p.adjust")
    p2 <- p2 + theme_opt
    
    p3 <- heatplot(gsea, 
                   showCategory = min(10, nrow(gsea)), 
                   foldChange = geneList)
    
    # Filter results for filtered terms and re-create plots
    filtered_result <- gsea@result %>% filter(ID %in% go_terms)
    gsea_filt <- gsea
    gsea_filt@result <- filtered_result
    
    p4 <- dotplot(gsea_filt, 
                  showCategory = min(10, nrow(filtered_result)))
    p4 <- p4 + theme_opt
    
    p5 <- ridgeplot(gsea_filt, 
                    showCategory = min(10, nrow(filtered_result)), 
                    fill = "p.adjust")
    p5 <- p5 + theme_opt
    
    p6 <- heatplot(gsea_filt, 
                   showCategory = min(10, nrow(filtered_result)), 
                   foldChange = geneList)
    p6 <- p6 + theme_opt
    
    p_all <- plot_grid(p1, p2, p3, p4, p5, p6, 
                       ncol = 1, 
                       labels = LETTERS[1:6], 
                       rel_heights = rep(1, 6))
    
    print("Plots generated")
    return(list(data = gsea@result, fig = p_all))
}

plot_expression_timecourse <- function(dds, sig_genes, 
                                       metaData, control, 
                                       padj_cutoff = 0.05,
                                       log2FC_cutoff = 1,
                                       use_rlog = TRUE) {
    
    # Extract normalized counts
    if (use_rlog) {
        norm_counts <- assay(rlog(dds, blind = FALSE))[sig_genes, ]
    } else {
        norm_counts <- counts(dds, normalized = TRUE)[sig_genes, ]
    }
    
    # Melt the data
    norm_counts_melt <- melt(norm_counts, 
                             varnames = c("gene", "sample"), 
                             value.name = "expression")
    
    # Merge with metadata
    norm_counts_melt <- merge(norm_counts_melt, 
                              metaData, 
                              by.x = "sample")
    norm_counts_melt$gene <- as.character(norm_counts_melt$gene)
    
    # Set factor levels for proper ordering
    treatments <- unique(metaData[["group"]][!grepl(control, 
                                                    metaData[["group"]])])
    group_levels <- c(control, treatments)
    norm_counts_melt$group <- factor(norm_counts_melt$group, 
                                     levels = group_levels)
    
    # Calculate mean and standard error
    norm_counts_summary <- norm_counts_melt %>%
        group_by(group, gene) %>%
        summarise(
            mean_expression = mean(expression),
            se_expression = sd(expression) / sqrt(n()),
            .groups = 'drop'
        ) %>% arrange(gene)

    # Create the plot
    p <- ggplot(norm_counts_summary, aes(x = group, 
                                         y = mean_expression, 
                                         color = gene, 
                                         group = gene)) + 
        geom_line() + 
        geom_point() + 
        geom_errorbar(aes(ymin = mean_expression - se_expression, 
                          ymax = mean_expression + se_expression), 
                      width = 0.2) + 
        facet_wrap(~ gene, scales = "free_y") + 
        theme_bw() + 
        labs(title = "Gene Expression Over Time", 
             x = "Condition", 
             y = "Normalized Expression") + 
        theme(axis.text.x = element_text(angle = 45, hjust = 1), 
              legend.position = "none")
    
    #p2 <- heatmap(norm_counts, group_levels)
    
    return(list(tc_plot = p1))#, hm = p2))
}


## Function to generate heatmap from rlog transformed data
rld_heatmap <- function(counts, groups) {
    annoDF <- data.frame(colData(rld)[, "group"])
    colnames(annoDF) <- "group"
    rownames(annoDF) <- colnames(rld)
    rld_gene_data <- assay(rld)[rownames(assay(rld)) %in% gene_list, ]
    rld_gene_data <- rld_gene_data[, order(colnames(rld_gene_data))]
    color_map <- c("#E69F00", "#56B4E9", "#009E73", 
                   "#F0E442", "#0072B2", "#D55E00", 
                   "#CC79A7", "#000000", "#999990")
    ann_colors = list(group = color_map[1:length(unique(metadata[["group"]]))])
    names(ann_colors$group) <- unique(metadata[["group"]])
    
    # Generate heatmap
    rld_hm <- pheatmap(
        rld_gene_data,
        scale = "row",
        show_rownames = TRUE,
        show_colnames = TRUE,
        treeheight_col = 2,
        treeheight_row = 0,
        cluster_rows = TRUE,
        cluster_cols = FALSE,
        annotation_col = annoDF,
        annotation_names_col = FALSE,
        annotation_colors = ann_colors,
        fontsize = 8,
        cellwidth = 10,
        cellheight = 2,
        border_color = NA
    )
    return(rld_hm)
}


## Use tryCatch for error handling
tryCatch({
    countData <- read.csv("read_counts.txt", header = TRUE, sep = "\t")
    metadata_file <- "sample_metadata.txt"
    metaData <- read.csv(sample_metadata.txt, header = TRUE, sep = "\t")
    go_terms <- readLines("immune_GO_numbers.txt")
    gene_list <- readLines("reactome_pathways/immune_system.tsv")
}, error = function(e) {
    stop("Error reading input files: ", e$message)
})


## Set variable values
control <- "HFF"
p_cutoff <- 0.05
fc_cutoff <- 1


## Set up DESeq2 dataset
dds <- suppressWarnings(DESeqDataSetFromMatrix(countData, 
                                               metaData, 
                                               as.formula("~group"), 
                                               tidy = TRUE))
dds[["group"]] <- relevel(dds[["group"]], ref = control)
dds <- DESeq(dds)


## Set up lists for excel workbooks and gene lists from DESeq2 results
wb <- createWorkbook()
results_list <- list()

for (treatment in unique(metaData[["group"]][!grepl(control, 
                                                    metaData[["group"]])])) {

    ## Perform DESeq2 comparison between control and treatment
    res <- results(dds, contrast = c("group", treatment, control))
    res <- as.data.frame(res)
    res <- res[order(res$log2FoldChange, decreasing = TRUE), ]
    res <- cbind(rownames(res), res)
    colnames(res)[colnames(res) == "rownames(res)"] <- "geneID"
    results_list[[treatment]] <- rownames(res)[which(abs(res$log2FoldChange) 
                                                     > fc_cutoff & res$padj 
                                                     < p_cutoff & rownames(res) 
                                                     %in% gene_list)]

    ## Perform GSEA and generate plots
    tryCatch({
        GSEA_results <- GSEA_analysis(res, 
                                      paste0(treatment, " vs ", control), 
                                      "ALL")
        
        if (!is.null(GSEA_results)) {
            pdf(file.path(dirname(metadata_file), 
                          paste0(treatment, "_vs_", control, 
                                 "_GSEA_plots.pdf")), 
                width = 10, 
                height = 30)
            print(GSEA_results$fig)
            dev.off()
            addWorksheet(wb, 
                         sheetName = paste0(treatment, "_vs_", control))
            writeData(wb, 
                      sheet = paste0(treatment, "_vs_", control), 
                      GSEA_results$data)
            print(paste("GSEA results for", treatment, "written to workbook"))
        } else {
            print(paste("No GSEA results for", treatment))
        }
    }, error = function(e) {
        print(paste("Error in GSEA analysis for", treatment, ":", e$message))
    })
}

saveWorkbook(wb, 
             file = file.path(dirname(metadata_file), 
                              paste0(control, "_control_GSEA_data.xlsx")), 
             overwrite = TRUE)
sig_genes <- unique(unlist(results_list))

## Perform rlog transformation and generate heatmap from genes in gene list
if (length(sig_genes) == 0) {
    warning("No genes pass the significance and fold change thresholds.")
    return(NULL)
} else {
    print(paste("Significant unique genes from gene list:", 
                paste(length(sig_genes), collapse = ", ")))
    tc_plot <- plot_expression_timecourse(dds, 
                                          sig_genes, 
                                          metaData, 
                                          control, 
                                          padj_cutoff = p_cutoff, 
                                          log2FC_cutoff = fc_cutoff, 
                                          use_rlog = TRUE)
    pdf(file.path(dirname(metadata_file), 
                  paste0(treatment, "_vs_", control, 
                         "_time_course_plot.pdf")), 
        width = 10, height = 8)
    print(tc_plot)
    dev.off()
}
