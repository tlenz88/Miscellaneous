# Load required libraries
library(DESeq2)
library(pheatmap)
library(gridExtra)
library(grid)
library(topGO)

# Load necessary data
countData <- read.csv("read_counts.txt", sep="\t")
metaData <- read.csv("sample_metadata.txt", sep="\t", row.names=1)
colnames(metaData) <- tolower(colnames(metaData))
GO_map <- readMappings("F:/organism_genome/Pfalciparum3D7_v66/PlasmoDB-66_Pfalciparum3D7_GO.map")
control <- "DMSO"
treatment <- "RAP"

# Set Gene_ID as row names for countData and remove it as a column
rownames(countData) <- countData$Gene_ID
countData$Gene_ID <- NULL

# Perform genome-wide differential expression analysis
dds <- DESeqDataSetFromMatrix(countData = countData, 
                              colData = metaData, 
                              design = ~group)
dds$group <- relevel(dds$group, ref = control)  # Assuming "WT" is the control group
dds <- DESeq(dds)
res <- results(dds, contrast = c("group", treatment, control))

# Get normalized counts
normalized_counts <- counts(dds, normalized=TRUE)

# Function for GO enrichment
GOenrichment <- function(DEGenes, gene_names, GO_map) {
    gene_list <- factor(as.integer(gene_names %in% DEGenes))
    names(gene_list) <- gene_names
    GO_data <- new("topGOdata",
                   ontology = "BP",
                   allGenes = gene_list,
                   annot = annFUN.gene2GO,
                   gene2GO = GO_map
    )
    fisher_test <- runTest(GO_data, statistic = "fisher")
    fisher_test_results <- GenTable(GO_data,
                                    weighted_Fisher = fisher_test,
                                    orderBy = "weighted_Fisher",
                                    ranksOf = "weighted_Fisher",
                                    topNodes = 50
    )
    fisher_terms = c(fisher_test_results$GO.ID)
    fisher_genes <- genesInTerm(GO_data, fisher_terms)
    fisher_pval <- as.vector(as.numeric(fisher_test_results$weighted_Fisher))
    fisher_test_results$logpval <- -log10(fisher_pval)
    return(fisher_test_results)
}

# Separate upregulated and downregulated genes
DEGenes_pos <- rownames(res)[which(res$padj < 0.05 & res$log2FoldChange > 0)]
DEGenes_neg <- rownames(res)[which(res$padj < 0.05 & res$log2FoldChange < 0)]

# Perform GO enrichment
gene_names <- rownames(res)
GO_results_pos <- GOenrichment(DEGenes_pos, gene_names, GO_map)
GO_results_neg <- GOenrichment(DEGenes_neg, gene_names, GO_map)

# Function to get genes for a specific GO term
getGenesForGOTerm <- function(GO_term, GO_map) {
    genes <- names(GO_map)[sapply(GO_map, function(x) GO_term %in% x)]
    return(genes)
}

# Function to create a heatmap
createHeatmap <- function(genes, normalized_counts, metaData, title, breaks) {
    genes <- genes[genes %in% rownames(normalized_counts)]
    
    if(length(genes) == 0) {
        warning(paste("No genes found for", title))
        return(NULL)
    }
    
    gene_counts <- normalized_counts[genes, , drop = FALSE]
    
    if(nrow(gene_counts) < 2 || ncol(gene_counts) < 2) {
        warning(paste("Not enough data to create heatmap for", title))
        return(NULL)
    }
    
    if(all(apply(gene_counts, 1, var) == 0)) {
        warning(paste("No variation in gene expression for", title))
        return(NULL)
    }
    
    mat <- log2(gene_counts + 1)  # log2 transform
    mat_scaled <- na.omit(t(scale(t(mat))))  # scale rows
    
    tryCatch({
        if (nrow(mat_scaled) > 30) {
            hm <- pheatmap(mat_scaled, 
                           main = title,
                           show_rownames = FALSE, 
                           show_colnames = FALSE, 
                           cluster_cols = FALSE,
                           annotation_col = metaData["group"],
                           breaks = breaks,
                           silent = TRUE)
        } else {
            hm <- pheatmap(mat_scaled, 
                           main = title,
                           show_rownames = TRUE, 
                           show_colnames = FALSE, 
                           cluster_cols = FALSE,
                           annotation_col = metaData["group"],
                           breaks = breaks,
                           silent = TRUE) 
        }    
        hm_grob <- grid.grabExpr(print(hm))
        return(hm_grob)
    }, error = function(e) {
        warning(paste("Error creating heatmap for", title, ":", e$message))
        return(NULL)
    })
}

# Get top 2 GO terms for upregulated and downregulated genes
top_GO_pos <- GO_results_pos$GO.ID[1:2]
top_GO_neg <- GO_results_neg$GO.ID[1:2]

# Create list of heatmaps
heatmap_list <- list()

# Function to add heatmap to list
add_heatmap <- function(genes, title, list_name, breaks) {
    hm <- createHeatmap(genes, normalized_counts, metaData, title, breaks)
    if (!is.null(hm)) {
        heatmap_list[[list_name]] <<- hm
    }
}

# Calculate global breaks for consistent color scale
all_genes <- unique(c(
    unlist(lapply(top_GO_pos, getGenesForGOTerm, GO_map)),
    unlist(lapply(top_GO_neg, getGenesForGOTerm, GO_map)),
    c("PF3D7_1463400", "PF3D7_1465200", "PF3D7_1342700", "PF3D7_0318200", 
      "PF3D7_0215700", "PF3D7_1364800", "PF3D7_1104700"),
    c("PF3D7_0309000", "PF3D7_1355500", "PF3D7_1135100", "PF3D7_0805600", 
      "PF3D7_1455100", "PF3D7_0802500", "PF3D7_0515900", "PF3D7_1322000", 
      "PF3D7_1469200", "PF3D7_1423300", "PF3D7_0810300", "PF3D7_0704800", 
      "PF3D7_1414400"),
    c("PF3D7_0321400", "PF3D7_1401800", "PF3D7_0217500", "PF3D7_0934800", 
      "PF3D7_0902100", "PF3D7_1450000", "PF3D7_0214600", "PF3D7_1223100", 
      "PF3D7_1124600", "PF3D7_1112100", "PF3D7_0203100", "PF3D7_1238900", 
      "PF3D7_1471400", "PF3D7_1148000", "PF3D7_1436600", "PF3D7_1414500", 
      "PF3D7_0311400", "PF3D7_0525900", "PF3D7_0928900", "PF3D7_1356800", 
      "PF3D7_0509800", "PF3D7_1113900")
))

all_genes <- all_genes[all_genes %in% rownames(normalized_counts)]
mat_all <- log2(normalized_counts[all_genes,] + 1)
mat_all_scaled <- na.omit(t(scale(t(mat_all))))
breaks <- seq(min(mat_all_scaled), max(mat_all_scaled), length.out = 100)

# Heatmaps for top 2 GO terms (upregulated)
for (i in 1:2) {
    genes <- getGenesForGOTerm(top_GO_pos[i], GO_map)
    add_heatmap(genes, paste0("Upregulated ", top_GO_pos[i]), paste0("GO_pos_", i), breaks)
}

# Heatmaps for top 2 GO terms (downregulated)
for (i in 1:2) {
    genes <- getGenesForGOTerm(top_GO_neg[i], GO_map)
    add_heatmap(genes, paste0("Downregulated ", top_GO_neg[i]), paste0("GO_neg_", i), breaks)
}

# DNA polymerases heatmap
dna_pol_genes <- c("PF3D7_1463400", "PF3D7_1465200", "PF3D7_1342700", "PF3D7_0318200", 
                   "PF3D7_0215700", "PF3D7_1364800", "PF3D7_1104700")
add_heatmap(dna_pol_genes, "DNA Polymerases", "DNA_polymerases", breaks)

# Phosphatases heatmap
phosphatase_genes <- c("PF3D7_0309000", "PF3D7_1355500", "PF3D7_1135100", "PF3D7_0805600", 
                       "PF3D7_1455100", "PF3D7_0802500", "PF3D7_0515900", "PF3D7_1322000", 
                       "PF3D7_1469200", "PF3D7_1423300", "PF3D7_0810300", "PF3D7_0704800", 
                       "PF3D7_1414400")
add_heatmap(phosphatase_genes, "Phosphatases", "Phosphatases", breaks)

# Kinases heatmap
kinase_genes <- c("PF3D7_0321400", "PF3D7_1401800", "PF3D7_0217500", "PF3D7_0934800", 
                  "PF3D7_0902100", "PF3D7_1450000", "PF3D7_0214600", "PF3D7_1223100", 
                  "PF3D7_1124600", "PF3D7_1112100", "PF3D7_0203100", "PF3D7_1238900", 
                  "PF3D7_1471400", "PF3D7_1148000", "PF3D7_1436600", "PF3D7_1414500", 
                  "PF3D7_0311400", "PF3D7_0525900", "PF3D7_0928900", "PF3D7_1356800", 
                  "PF3D7_0509800", "PF3D7_1113900")
add_heatmap(kinase_genes, "Kinases", "Kinases", breaks)

# Combine all heatmaps into a single figure
if (length(heatmap_list) > 0) {
    pdf("combined_heatmaps.pdf", width = 20, height = 30)
    do.call(grid.arrange, c(heatmap_list, ncol = 2))
    dev.off()
} else {
    warning("No heatmaps were created. Check your data and gene lists.")
}