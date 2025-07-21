Quality_assessment <- function(metadata){
    library(BiocManager)
    library(ChIPQC)
    
    # Read peaksets
    samples <- read.csv(metadata, sep='\t')
    df <- dba(sampleSheet=samples)
    df <- dba.count(df, bUseSummarizeOverlaps=TRUE)
    
    # Normalize data
    df <- dba.normalize(df)
    norm <- dba.normalize(df, bRetrieve=TRUE)
    normlibs <- cbind(FullLibSize=norm$lib.sizes, NormFacs=norm$norm.factors,
                      NormLibSize=round(norm$lib.sizes/norm$norm.factors))
    rownames(normlibs) <- info$ID
    
    # Establish model design and contrast
    df <- dba.contrast(df, reorderMeta=list(Condition="control"))
    
    # Perform differential analysis
    df <- dba.analyze(df)
    dba.show(df, bContrasts=TRUE)
    
    # Retrieve differentially bound sites
    df.DB <- dba.report(df)
    
    return(df)
}

Quality_plot <- function(df, pdf_output){
    pdf(pdf_output)
    
    # Correlation Heatmap
    corrvals <- dba.plotHeatmap(df)
    
    # PCA plot
    dba.plotPCA(df, contrast=1, method=DBA_DESEQ2, attributres=DBA_FACTOR, label=DBA_CONDITION)
    
    # Venn diagram
    dba.plotVenn(df, contrast=1, method=DBA_DESEQ2)
    
    # MA plot
    dba.plotMA(df, dotSize=1, bSmooth=FALSE)
    dba.plotMA(df, dotSize=1, bSmooth=FALSE, bXY=TRUE)
    
    # Boxplot
    dba.plotBox(df)
    
    # Volcano plot
    dba.plotVolcano(df)
    
    # Heatmaps
    hmap <- colorRampPalette(c("blue", "white", "red"))(n = 15)
    readscores <- dba.plotHeatmap(df, contrast=1, correlations=FALSE, scale="row", colScheme=hmap)
    
    dev.off()
}
