required_packages <- c("DESeq2", "ggplot2", "pheatmap", "RColorBrewer", "GenomicRanges",
                       "apeglm", "dplyr", "EnhancedVolcano", "readr", "rlang", "fgsea",
                       "rlang")

missing_packages <- required_packages[!sapply(required_packages, requireNamespace, quietly = TRUE)]

if (length(missing_packages) > 0) {
    cat("Error: The following required packages are missing:\n")
    cat(paste(missing_packages, collapse = ", "), "\n")
    cat("Please install them before running this script.\n")
    quit(status = 1)
}

# Load packages
suppressPackageStartupMessages({
    invisible(lapply(required_packages, library, character.only = TRUE))
})

# Define path to the Reactome pathway files
pathway_dir <- "F:/toxo_project/RNAseq/Hsapien_output/reactome_pathways/"
pathway_files <- list.files(pathway_dir, pattern = "*.tsv", full.names = TRUE)
pathway_names <- gsub("\\.tsv$", "", basename(pathway_files))

# Read pathways into a named list of data frames
gene_pathway_dataframes <- lapply(pathway_files, function(file) read.csv(file, sep="\t", header=FALSE))
names(gene_pathway_dataframes) <- pathway_names

# Read files and define input variables
read_counts <- "F:/toxo_project/RNAseq/Hsapien_output/DEGenes_24ACD/read_counts.txt"
metadata <- "F:/toxo_project/RNAseq/Hsapien_output/DEGenes_24ACD/sample_metadata.txt"
countData <- read.csv(read_counts, header = TRUE, sep = "\t")
metaData <- read.csv(metadata, header = TRUE, sep = "\t")
colnames(metaData) <- tolower(colnames(metaData))
group_column <- "group"
qval <- 0.05
control <- "HFF"
num_top_genes <- 2000
key_pathway <- "immune_system"

# Create DESeq object and rlog transformation
dds <- DESeqDataSetFromMatrix(countData, metaData, as.formula(paste0("~", group_column)), tidy = TRUE)
dds[[group_column]] <- relevel(dds[[group_column]], ref = control)
dds <- DESeq(dds)
rld <- rlog(dds, blind = FALSE)

wong_colormap <- c("#E69F00", "#56B4E9", "#CC79A7", "#F0E442", "#009E73",  "#0072B2", "#D55E00")

for (topVarGenes in c(1000, 2000, 3000, 4000, 5000)) {
    ranked_genes <- rowVars(assay(rld))
    names(ranked_genes) <- rownames(rld)
    ranked_genes <- sort(ranked_genes, decreasing = TRUE)[1:topVarGenes]
    
    # Create dataframe for plotting
    topVarGenes_index <- head(order(rowVars(assay(rld)), decreasing = TRUE), topVarGenes)
    sigMat <- assay(rld)[topVarGenes_index, ]
    sigMat <- sigMat[, order(colnames(sigMat))]
    
    # Create a pathways list from your gene_pathway_dataframes
    pathways_list <- lapply(gene_pathway_dataframes, function(df) df[, 1])
    
    # Run GSEA with fgsea
    invisible(capture.output({
        gsea_results <- fgsea(pathways = pathways_list, stats = ranked_genes, minSize = 10, maxSize = 3000, scoreType = "pos", nproc = 1)
    }))
    
    # Filter significant pathways based on an adjusted p-value threshold
    significant_pathways <- data.frame(gsea_results[gsea_results$padj < 0.05, ])
    
    # Generate colormap for plotting using visually distinct color scale
    pathway_colormap <- wong_colormap[1:nrow(significant_pathways)]
    names(pathway_colormap) <- significant_pathways$pathway
    
    cellheight = 0.39 * 1000 / topVarGenes
    
    # Initialize rowAnno with FALSE and annotated_pathway with NA
    rowAnno <- data.frame(matrix(FALSE, nrow = nrow(sigMat), ncol = nrow(significant_pathways)))
    colnames(rowAnno) <- significant_pathways$pathway
    rownames(rowAnno) <- rownames(sigMat)
    
    sigMat_annotated <- as.data.frame(sigMat)
    sigMat_annotated$annotated_pathway <- rep(NA_character_, nrow(sigMat_annotated))
    sigMat_annotated$pathway_color <- rep(NA_character_, nrow(sigMat_annotated))
    
    # Annotate each gene based on pathway membership in significant pathways
    for (i in seq_len(nrow(significant_pathways))) {
        pathway_name <- significant_pathways$pathway[i]
        pathway_genes <- gene_pathway_dataframes[[pathway_name]][, 1]
        
        # Mark genes in rowAnno based on pathway membership
        rowAnno[, pathway_name] <- rownames(sigMat) %in% pathway_genes
    }
    
    # Assign the pathway with the lowest padj to annotated_pathway for each gene
    for (gene in rownames(sigMat_annotated)) {
        gene_pathways <- significant_pathways$pathway[unlist(rowAnno[gene, ])]
        
        if (length(gene_pathways) > 0) {
            padj_values <- significant_pathways[significant_pathways$pathway %in% gene_pathways, "padj"]
            padj_values <- unlist(padj_values)
            
            # Find the pathway with the lowest padj for this gene
            padj_pathways <- significant_pathways[significant_pathways$pathway %in% gene_pathways, ]
            
            # Check if there's another pathway with the same padj as key_pathway
            if (key_pathway %in% padj_pathways$pathway) {
                min_padj_pathway <- key_pathway
            } else {
                min_padj_pathway <- padj_pathways[order(padj_pathways$padj), "pathway"][1]
            }
            
            sigMat_annotated[gene, "annotated_pathway"] <- min_padj_pathway
            sigMat_annotated[gene, "pathway_color"] <- pathway_colormap[min_padj_pathway]
        }
    }

    if (exists("key_pathway") && key_pathway %in% significant_pathways$pathway) {
        break
    }
}


# Create a data frame for row annotations, with the pathway color for each gene
rowAnno <- data.frame(annotated_pathway = sigMat_annotated$annotated_pathway)

# Make sure the row names of row_anno match the row names of sigMat_annotated
rownames(rowAnno) <- rownames(sigMat_annotated)

# Prepare the expression data without the annotated columns
expression_data <- sigMat_annotated[, 1:(ncol(sigMat_annotated) - 2)]

colAnno <- data.frame(colData(rld)[, group_column])
colnames(colAnno) <- "group"
rownames(colAnno) <- colnames(rld)
group_colormap <- c("#5A5A5A", "#000000", "#D9D9D9")
names(group_colormap) <- unique(colAnno$group)

colorAnno <- list(group = group_colormap, annotated_pathway = pathway_colormap)

pdf("F:/toxo_project/RNAseq/Hsapien_output/test_reactome_heatmaps.pdf")
pheatmap(
    expression_data,
    main = paste("Top", topVarGenes, "most variable genes among all samples"),
    scale = "row",
    show_rownames = FALSE,
    show_colnames = TRUE,
    treeheight_col = 2,
    treeheight_row = 0,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    annotation_col = colAnno,
    annotation_row = rowAnno,
    annotation_names_col = FALSE,
    annotation_colors = colorAnno,
    fontsize = 8,
    cellwidth = 18,
    cellheight = cellheight,
    border_color = NA
)
dev.off()




VolcanoPlot <- function(df, pathway_file, treatment, control, min_log2FC, max_log2FC, min_padj, max_padj) {
    pathway <- read.csv(paste0("F:/toxo_project/RNAseq/Hsapien_output/reactome_pathways/", pathway_file, ".tsv"), header=FALSE)
    pathway_df <- df[which(rownames(df) %in% pathway[,1]), ]
    pathway_df_sig <- pathway_df[pathway_df$padj < 0.05, ]
    pathway_df_sig_up <- pathway_df_sig[pathway_df_sig$log2FoldChange > 0, ]
    pathway_df_sig_down <- pathway_df_sig[pathway_df_sig$log2FoldChange < 0, ]
    df$color <- ifelse(
        abs(df$log2FoldChange) >= 1 & df$padj <= qval,
        ifelse(rownames(df) %in% pathway[,1], "#F3766E","#999999"), "#000000")
    p <- ggplot(df, aes(x = log2FoldChange, y = -log10(padj))) +
        geom_point(aes(color = color), size = 1, alpha = 0.5) +
        scale_color_identity() +
        labs(
            title = paste(treatment, "vs", control, "-", pathway_file),
            subtitle = paste(pathway_file, "genes"),
            x = "Log2 Fold Change",
            y = "-log10(padj)",
            caption = paste(nrow(pathway_df_sig), "genes in pathway with padj <=", qval, ",", nrow(pathway_df_sig_up), "genes upregulated (log2FC >= 1),", nrow(pathway_df_sig_down), "genes downregulated (log2FC <= -1)")
        ) +
        xlim(c(min_log2FC*1.1, max_log2FC*1.1)) +
        ylim(c(min_padj, max_padj)) +
        theme_minimal() +
        theme(
            axis.title = element_text(size = 8),
            axis.text = element_text(size = 8),
            legend.position = "none"
        )
    if (pathway_file == "immune_system") {
        p <- p +
            geom_point(data = df %>% filter(color == "#56B4E9"), aes(x = log2FoldChange, y = -log10(padj)), color = "#56B4E9", size = 2)
    } else if (pathway_file == "cellular_response_to_stimuli") {
        p <- p +
            geom_point(data = df %>% filter(color == "#E69F00"), aes(x = log2FoldChange, y = -log10(padj)), color = "#E69F00", size = 2)
    } else if (pathway_file == "signal_transduction") {
        p <- p +
            geom_point(data = df %>% filter(color == "#CC79A7"), aes(x = log2FoldChange, y = -log10(padj)), color = "#CC79A7", size = 2)
    }
    print(p)
}


VolcanoPlot <- function(df, pathway_file, treatment, control, min_log2FC, max_log2FC, min_padj, max_padj, qval = 0.05) {
    # Load pathway genes
    pathway <- read.csv(paste0("F:/toxo_project/RNAseq/Hsapien_output/reactome_pathways/", pathway_file, ".tsv"), header=FALSE)
    
    # Subset df to include only genes in the pathway
    pathway_df <- df[which(rownames(df) %in% pathway[,1]), ]

    # Define significant upregulated and downregulated genes in the pathway
    pathway_df_sig <- pathway_df[pathway_df$padj < qval, ]
    pathway_df_sig_up <- pathway_df_sig[pathway_df_sig$log2FoldChange > 0, ]
    pathway_df_sig_down <- pathway_df_sig[pathway_df_sig$log2FoldChange < 0, ]
    
    # Assign colors based on significance within the pathway genes
    if (pathway_file == "immune_system") {
        pathway_df$color <- ifelse(abs(pathway_df$log2FoldChange) >= 1 & pathway_df$padj <= qval, "#56B4E9", "#000000")
    } else if (pathway_file == "cellular_response_to_stimuli") {
        pathway_df$color <- ifelse(abs(pathway_df$log2FoldChange) >= 1 & pathway_df$padj <= qval, "#E69F00", "#000000")
    } else if (pathway_file == "signal_transduction") {
        pathway_df$color <- ifelse(abs(pathway_df$log2FoldChange) >= 1 & pathway_df$padj <= qval, "#CC79A7", "#000000")
    }
    
    # Create volcano plot with only pathway genes
    p <- ggplot(pathway_df, aes(x = log2FoldChange, y = -log10(padj))) +
        geom_point(aes(color = color), size = 1, alpha = 0.5) +
        scale_color_identity() +
        labs(
            title = paste(treatment, "vs", control, "-", pathway_file),
            subtitle = paste(pathway_file, "genes"),
            x = "Log2 Fold Change",
            y = "-log10(padj)",
            caption = paste(nrow(pathway_df_sig), "genes in pathway with padj <=", qval, ",", 
                            nrow(pathway_df_sig_up), "genes upregulated (log2FC >= 1),", 
                            nrow(pathway_df_sig_down), "genes downregulated (log2FC <= -1)")
        ) +
        xlim(c(min_log2FC * 1.1, max_log2FC * 1.1)) +
        ylim(c(min_padj, max_padj)) +
        theme_minimal() +
        theme(
            axis.title = element_text(size = 8),
            axis.text = element_text(size = 8),
            legend.position = "none"
        )
    print(p)
}



pdf(paste0("F:/toxo_project/RNAseq/Hsapien_output/pathway_volcano_plots.pdf"))
for (control in c("HFF", "ME49RFP_6hpi")) {
    dds[[group_column]] <- relevel(dds[[group_column]], ref = control)
    dds <- DESeq(dds)
    for (treatment in c("ME49RFP_6hpi", "ME49RFP_24hpi")) {
        if (control == treatment) {
            next
        } else if ((control == "HFF") & (treatment == "ME49RFP_24hpi")) {
            next
        }
        resLFCapeglm <- lfcShrink(
            dds,
            coef = paste0("group_", treatment, "_vs_", control),
            type = "apeglm"
        )
        resLFCapeglm <- resLFCapeglm[!is.na(resLFCapeglm$log2FoldChange) & !is.na(resLFCapeglm$padj), ]
        #min_log2FC <- min(resLFCapeglm$log2FoldChange, na.rm = TRUE)
        min_log2FC <- -6.1
        #max_log2FC <- max(resLFCapeglm$log2FoldChange, na.rm = TRUE)
        max_log2FC <- 11
        #min_padj <- -log10(max(resLFCapeglm$padj, na.rm = TRUE))
        min_padj <- 0
        #max_padj <- -log10(min(resLFCapeglm$padj, na.rm = TRUE))
        max_padj <- 30
        df <- as.data.frame(resLFCapeglm)
        df$gene <- rownames(df)
        for (pathway_file in c("immune_system", "cellular_response_to_stimuli", "signal_transduction")) {
            VolcanoPlot(df, pathway_file, treatment, control, min_log2FC, max_log2FC, min_padj, max_padj)
        }
    }
}
dev.off()