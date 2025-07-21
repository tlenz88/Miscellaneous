df <- read.csv(paste0("F:/Deitsch_collab/HiC/matrices/merged/D2/test.txt"), header = FALSE, sep = "\t")
colnames(df) <- c("i","j","q")

mat <- matrix(0.05, 2337,2337)

for(i in 1:nrow(df)) {
  mat[df[i,1],df[i,2]] <- df[i,3]
  mat[df[i,2],df[i,1]] <- df[i,3]
}

mat <- apply(t(mat),2,rev)

pheatmap(mat, 
         main = "D2 10kb significant\ninterchromosomal interactions",
         show_colnames = FALSE, 
         show_rownames = FALSE, 
         cluster_cols = FALSE, 
         cluster_rows = FALSE, 
         colorRampPalette(rev(brewer.pal(9, "Blues")))(255), 
         cellwidth = 360 / ncol(mat), 
         cellheight = 360 / ncol(mat),
         legend_breaks = c(0, 0.01, 0.02, 0.03, 0.04, 0.05))
