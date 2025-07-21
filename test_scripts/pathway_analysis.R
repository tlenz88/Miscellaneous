#!/usr/bin/env Rscript

## Created: September 5, 2024
## Updated: September 5, 2024
## Author(s): Todd Lenz, tlenz001@ucr.edu

## Performs pathway analysis.

args <- commandArgs(trailingOnly = TRUE)

ShowHelp <- function() {
    cat("differential_expression.R --help\n")
    cat("Usage: differential_expression.R -m METADATA -c CONTROL -q QVALUE [OPTIONS]\n")
    cat("-------------------------------------------------------------------\n")
    cat("  Input arguments:\n")
    cat("    -m|--metadata METADATA : Metadata file.\n")
    cat("    -c|--control CONTROL   : String indicating control group.\n")
    cat("  Optional arguments:\n")
    cat("    -g|--group_col GROUP   : Group column in metadata (default: 'group').\n")
    cat("    -q|--qvalue QVALUE     : Q-value cutoff (default: 0.05).\n")
    cat("    -h|--help HELP         : Show help message.\n")
    cat("-------------------------------------------------------------------\n")
}

if ("--help" %in% args || "-h" %in% commandArgs() || length(commandArgs(trailingOnly = TRUE)) == 0) {
    ShowHelp()
    quit()
}

###########################################
## Load packages and install if missing. ##
###########################################
required_packages <- c("org.Hs.eg.db", "clusterProfiler", "ReactomePA", "ggplot2", "openxlsx")

# Function to install CRAN packages
install_cran_packages <- function(packages) {
    # Identify missing packages
    missing_packages <- packages[!(packages %in% installed.packages()[, "Package"])]
    
    # Install missing CRAN packages
    if (length(missing_packages) > 0) {
        message("Installing missing CRAN packages: ", paste(missing_packages, collapse = ", "))
        install.packages(missing_packages, repos = "https://cloud.r-project.org/")
    }
    
    # Load packages with suppressed messages
    invisible(lapply(packages, function(pkg) {
        suppressMessages(require(pkg, character.only = TRUE))
    }))
}

# Install CRAN packages
install_cran_packages(required_packages)

# Install Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager", repos = "https://cloud.r-project.org/")
}

# Function to install missing Bioconductor packages
install_bioc_packages <- function(packages) {
    # Identify missing Bioconductor packages
    missing_packages <- packages[!(packages %in% installed.packages()[, "Package"])]
    
    # Install missing Bioconductor packages
    if (length(missing_packages) > 0) {
        message("Installing missing Bioconductor packages: ", paste(missing_packages, collapse = ", "))
        BiocManager::install(missing_packages, update = FALSE, ask = FALSE)
    }
}

# Install missing Bioconductor packages
install_bioc_packages(required_packages)

# Ensure all required packages are loaded with suppressed messages
invisible(lapply(required_packages, function(pkg) {
    suppressMessages(requireNamespace(pkg, quietly = TRUE))
}))

############################
## Check input arguments. ##
############################
parse_args <- function(args) {
    parsed <- list()
    i <- 1
    while (i <= length(args)) {
        if (args[i] %in% c("--metadata", "-m")) {
            parsed$metadata <- args[i + 1]
            i <- i + 2
        } else if (args[i] %in% c("--control", "-c")) {
            parsed$control <- args[i + 1]
            i <- i + 2
        } else if (args[i] == "--group_column") {
            parsed$group_column <- args[i + 1]
            i <- i + 2
        } else if (args[i] %in% c("--qval", "-q")) {
            parsed$qval <- as.numeric(args[i + 1])
            i <- i + 2
        } else {
            cat("Error: Unknown argument", args[i], "\n")
            quit(status = 1)
        }
    }
    return(parsed)
}

parsed_args <- parse_args(args)

control <- parsed_args$control
if (is.null(parsed_args$group_column)) parsed_args$group_column <- "group"
if (is.null(parsed_args$qval)) parsed_args$qval <- 0.05

###########################
## Load sample metadata. ##
###########################
safe_read_csv <- function(file_path, ...) {
    tryCatch(
        {
            read.csv(file_path, ...)
        },
        error = function(e) {
            cat("Error reading file:", file_path, "\n")
            cat("Error message:", e$message, "\n")
            quit(status = 1)
        }
    )
}

metaData <- safe_read_csv(parsed_args$metadata, header = TRUE, sep = "\t")
colnames(metaData) <- tolower(colnames(metaData))

###############################
## Perform pathway analysis. ##
###############################
PathwayAnalysis <- function(deg, title, treatment) {
    if (nrow(deg) == 0) {
        message("No genes passed the filtering criteria.")
        return(NULL)
    }
    entrezIDs <- mapIds(org.Hs.eg.db, keys = deg$geneID, column = "ENTREZID", keytype = "SYMBOL")
    deg$entrezID <- entrezIDs
    deg <- na.omit(deg)
    deg_df <- as.data.frame(deg)
    missingIDs <- (1 - (length(na.omit(entrezIDs)) / length(entrezIDs))) * 100
    sprintf("%.2f%% of genes could not be mapped to ENTREZ IDs.", missingIDs)
    if (nrow(deg) == 0) {
        message("No Entrez IDs found for the given genes.")
        return(NULL)
    }
    reactome_enrich <- enrichPathway(gene = deg$entrezID,
                                     organism = "human",
                                     pvalueCutoff = 0.05,
                                     pAdjustMethod = "BH",
                                     qvalueCutoff = 0.2,
                                     readable = TRUE)
    if (is.null(reactome_enrich) || nrow(reactome_enrich) == 0) {
        message("No enriched pathways found.")
        return(NULL)
    }

    reactome_enrich_df <- as.data.frame(reactome_enrich)

    p <- dotplot(reactome_enrich, showCategory = 15, title = title, font.size = 8)
    p <- p + 
        theme_bw() + 
        theme(axis.text.y = element_text(size = 8, vjust = 0.5),
              axis.text.x = element_text(size = 8),
              axis.title.y = element_text(size=8),
              axis.title.x = element_text(size=8))
    print(p)
    return(list(df = reactome_enrich_df, deg_df = deg_df, plot = p))
}

wb <- createWorkbook()
data_frames_list <- list()

pdf(file.path(dirname(parsed_args$metadata), paste0(control, "_control_pathway_enrichment.pdf")))
for (treatment in unique(metaData[[parsed_args$group_column]][!grepl(control, metaData[[parsed_args$group_column]])])) {
    DEGenes <- file.path(dirname(parsed_args$metadata), paste0(treatment, "_vs_", control, "_DEgenes.txt"))
    deg <- read.csv(DEGenes, sep='\t')
    deg_up <- subset(deg, padj < parsed_args$qval & log2FoldChange > 1)
    deg_down <- subset(deg, padj < parsed_args$qval & log2FoldChange < -1)
    title_up <- paste0("Pathways upregulated in ", treatment, " vs ", control)
    title_down <- paste0("Pathways downregulated in ", treatment, " vs ", control)
    result_up <- PathwayAnalysis(deg_up, title_up, treatment)
    result_down <- PathwayAnalysis(deg_down, title_down, treatment)
    if (!is.null(result_up)) {
        data_frames_list[[paste0(treatment, "_upreg")]] <- result_up$deg_df
        data_frames_list[[paste0(treatment, "_pathway_upreg")]] <- result_up$df
    }
    if (!is.null(result_down)) {
        data_frames_list[[paste0(treatment, "_downreg")]] <- result_down$deg_df
        data_frames_list[[paste0(treatment, "_pathway_downreg")]] <- result_down$df
    }
}

for (sheet_name in names(data_frames_list)) {
    addWorksheet(wb, sheetName = sheet_name)
    writeData(wb, sheet = sheet_name, data_frames_list[[sheet_name]])
}

saveWorkbook(wb, file = file.path(dirname(parsed_args$metadata), paste0(control, "_control_pathway_enrichment.xlsx")), overwrite = TRUE)
