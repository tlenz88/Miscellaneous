library(dplyr)

chess_bed <- function(mat, bed) {
  colnames(mat) <- c('i', 'j', 'int')
  colnames(bed) <- c('chr1', 'start1', 'end1', 'id')
  new_mat <- left_join(mat, bed, by = c('i' = 'id'))
  colnames(bed) <- c('chr2', 'start2', 'end2', 'id')
  new_mat <- left_join(new_mat, bed, by = c('j' = 'id'))
  new_mat$id <- seq.int(1:nrow(new_mat))
  new_mat$misc <- "."
  new_mat$strand1 <- "+"
  new_mat$strand2 <- "+"
  new_mat <- new_mat[,4:13]
  return(new_mat)
}

mat <- read.table(paste0("F:/ap2_hic/matrices/AP2_16A_CTRL/AP2_16A_CTRL_10000_iced_noAPIorMIT.matrix"))
bed <- read.table(paste0("F:/ap2_hic/matrices/test/Pfalciparum3D7_200kbWindow_10kbStep.txt"))
new_mat <- chess_bed(mat, bed)


write.table(new_mat, "F:/ap2_hic/matrices/test/Pfalciparum3D7_10000_noAPIorMIT.bed", sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)