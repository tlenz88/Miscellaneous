library(ggplot2)
library(EnhancedVolcano)
library(pheatmap)
library(dplyr)
library(stringr)
library(RColorBrewer)

#################################################################################################################################
############################### INCREASE AND DESCREASE IN VAR GENE INTERACTIONS BETWEEN WT AND KO ###############################
#################################################################################################################################

var <- read.csv('F:/organism_genome/Pfalciparum3D7/var_genes.gff', header = FALSE, sep = "\t")[,c(1,4,5,9)]
colnames(var) <- c("chr", "bin1", "end", "gene")
for(i in 1:nrow(var)) {
  var[i, 2] <- ceiling(var[i, 2] / 10000)
  var[i, 3] <- ceiling(var[i, 3] / 10000)
}
var$gene <- str_extract_all(var$gene, "(?<==).*(?=;)")
var <- var %>% mutate_at("gene", as.character)

sizes <- read.csv('F:/organism_genome/Pfalciparum3D7/PlasmoDB-58_Pfalciparum3D7_Genome.chrom.sizes', header = FALSE, sep = "\t")[1:14,]
starts <- data.frame(matrix(NA, nrow = nrow(sizes), ncol = ncol(sizes)))
colnames(starts) <- c("chr", "start_bin")
for (i in 1:nrow(sizes)) {
  starts[i, 1] <- sizes[i, 1]
  starts[1, 2] <- 0
  if (i > nrow(sizes) - 1) {
    sizes[i, 2] <- ceiling(sizes[i, 2] / 10000)
    break
  } else {
    starts[i+1, 2] <- ceiling(sizes[i, 2] / 10000) + starts[i, 2]
    sizes[i, 2] <- ceiling(sizes[i, 2] / 10000)
    
  }
}

var <- left_join(starts, var, by = "chr")
var["bin1"] <- var["bin1"] + var["start_bin"]
var["end"] <- var["end"] + var["start_bin"]
var <- subset(var, select = -start_bin)
var <- var[!is.na(var$bin1),]

ko <- c("CSA")
for (i in ko) {
    df = read.csv(paste0("F:/Ribacke_collab/HiC/output_files/",i,"/",i,"_10000_iced_cpm.matrix"), header = FALSE, sep = "\t")
    colnames(df) = c("bin1", "bin2", "ko_int")
    assign(paste0(i), df)
}

CSA$ko_int <- CSA$ko_int / sum(CSA$ko_int) * 1000000

wt <- c("NF54")
for (i in wt) {
    df = read.csv(paste0("F:/Ribacke_collab/HiC/output_files/",i,"/",i,"_10000_iced.matrix"), header = FALSE, sep = "\t")
    colnames(df) = c("bin1", "bin2", "wt_int")
    assign(paste0(i), df)
}

NF54$wt_int <- NF54$wt_int / sum(NF54$wt_int) * 1000000

for (i in c("CSA", "NF54")) {
  df <- get(i)
  v <- left_join(var, df, by = "bin1")
  assign(paste0(i, "_var"), v)
}

for (i in c("CSA_var")) {
    df <- get(i)
    df_mean <- data.frame(tapply(df$ko_int, factor(df$gene), mean))
    df_mean$gene <- row.names(df_mean)
    colnames(df_mean) <- c("ko_mean", "gene")
    df_mean <- df_mean[, c("gene", "ko_mean")]
    assign(paste0(i, "_mean"), df_mean)
}

for (i in c("NF54_var")) {
    df <- get(i)
    df_mean <- data.frame(tapply(df$wt_int, factor(df$gene), mean))
    df_mean$gene <- row.names(df_mean)
    colnames(df_mean) <- c("wt_mean", "gene")
    df_mean <- df_mean[, c("gene", "wt_mean")]
    assign(paste0(i, "_mean"), df_mean)
}

NF54_CSA_var_mean <- left_join(NF54_var_mean, CSA_var_mean, by = "gene")

NF54_CSA_var <- left_join(NF54_var, CSA_var, by=c("chr", "bin1", "end", "gene", "bin2"))

for (i in c("NF54_CSA_var")) {
  df = get(i)
  df[is.na(df)] = 0
  df$change <- df$ko_int > df$wt_int
  df["change"][df["change"] == "TRUE"] <- "up"
  df["change"][df["change"] == "FALSE"] <- "down"
  df$change <- factor(df$change, levels = c("up", "down"))
  assign(paste0(i, "_int"), df)
}

twoSD <- max(2*sd(NF54_CSA_var_int$wt_int),
             2*sd(NF54_CSA_var_int$ko_int))

ggplot(NF54_CSA_var_int, aes(x = wt_int, y = ko_int, color = change)) + 
    geom_point(size = 2) + 
    ggtitle("CSA vs NF54") + 
    labs(y = "KO interactions", x = "WT interactions") + 
    scale_x_continuous(expand = c(0, 0), limits = c(-.02*twoSD, twoSD)) + 
    scale_y_continuous(expand = c(0, 0), limits = c(-.02*twoSD, twoSD)) + 
    geom_abline(intercept = 0, slope = 1, size = 1,
                caption = paste0("Increased var interactions: ", length(which(NF54_CSA_var_int$change=="up")), 
                                 "  Decreased var interactions: ", length(which(NF54_CSA_var_int$change=="down"))))

pdf("F:/Deitsch_collab/HiC/figures/scatter_plots/log2change_var.pdf")

dev.off()

#################################################################################################################################

var_ratio <- data.frame(c(length(which(A3A_B10A_var_log2change$change=="up")) / 
                            (length(which(A3A_B10A_var_log2change$change=="down")) + 
                               length(which(A3A_B10A_var_log2change$change=="up"))), 
                          length(which(A3A_D2A_var_log2change$change=="up")) / 
                            (length(which(A3A_D2A_var_log2change$change=="down")) + 
                               length(which(A3A_D2A_var_log2change$change=="up"))), 
                          length(which(A3A_G91A_var_log2change$change=="up")) / 
                            (length(which(A3A_G91A_var_log2change$change=="down")) + 
                               length(which(A3A_G91A_var_log2change$change=="up"))), 
                          length(which(A3B_B10B_var_log2change$change=="up")) / 
                            (length(which(A3B_B10B_var_log2change$change=="down")) + 
                               length(which(A3B_B10B_var_log2change$change=="up"))), 
                          length(which(A3B_D2C_var_log2change$change=="up")) / 
                            (length(which(A3B_D2C_var_log2change$change=="down")) + 
                               length(which(A3B_D2C_var_log2change$change=="up"))), 
                          length(which(A3B_G91B_var_log2change$change=="up")) / 
                            (length(which(A3B_G91B_var_log2change$change=="down")) + 
                               length(which(A3B_G91B_var_log2change$change=="up"))), 
                          length(which(A3A_B10A_var_log2change$change=="down")) / 
                            (length(which(A3A_B10A_var_log2change$change=="up")) + 
                               length(which(A3A_B10A_var_log2change$change=="down"))), 
                          length(which(A3A_D2A_var_log2change$change=="down")) / 
                            (length(which(A3A_D2A_var_log2change$change=="up")) + 
                               length(which(A3A_D2A_var_log2change$change=="down"))), 
                          length(which(A3A_G91A_var_log2change$change=="down")) / 
                            (length(which(A3A_G91A_var_log2change$change=="up")) + 
                               length(which(A3A_G91A_var_log2change$change=="down"))), 
                          length(which(A3B_B10B_var_log2change$change=="down")) / 
                            (length(which(A3B_B10B_var_log2change$change=="up")) + 
                               length(which(A3B_B10B_var_log2change$change=="down"))), 
                          length(which(A3B_D2C_var_log2change$change=="down")) / 
                            (length(which(A3B_D2C_var_log2change$change=="up")) + 
                               length(which(A3B_D2C_var_log2change$change=="down"))), 
                          length(which(A3B_G91B_var_log2change$change=="down")) / 
                            (length(which(A3B_G91B_var_log2change$change=="up")) + 
                               length(which(A3B_G91B_var_log2change$change=="down")))))

colnames(var_ratio) <- "ratio"
var_ratio$sample <- c("KO", "KO", "KO", "KO", "KO", "KO", "WT", "WT", "WT", "WT", "WT", "WT", "WT", "WT")
var_ratio$samples <- c("A3A vs B10A", "A3A vs D2A", "A3A vs G91A",
                       "A3B vs B10B", "A3B vs D2A", "A3B vs G91B",
                       "A3A vs B10A", "A3A vs D2A", "A3A vs G91A",
                       "A3B vs B10B", "A3B vs D2A", "A3B vs G91B")

pdf("D:/Deitsch_collab/figures/barplots/var_WT_vs_KO.pdf")
ggplot(var_ratio, aes(y = ratio, x = samples, fill = sample)) + 
  geom_bar(position = "stack", stat = "identity") + 
  geom_hline(yintercept = 0, size = 1) + 
  theme(axis.title.x = element_blank(), 
        axis.title.y = element_text(size = 12, face = "bold"), 
        axis.text.x = element_text(face = "bold", vjust = 10, size = 10), 
        axis.ticks.x = element_blank(), 
        axis.text.y = element_text(face = "bold", size = 10), 
        axis.ticks.y = element_line(size = 1), 
        axis.ticks.length.y = unit(.25, "cm"), 
        panel.background = element_blank(), 
        axis.line.y.left = element_line(size = 1, color = "black"))
dev.off()

#################################################################################################################################

ggplot(var_df, aes(y = change, x = samples, fill = ), geom_bar(stat = "identity", color = "blue"))

pdf("D:/Deitsch_collab/figures/scatter_plots/log2change_var2.pdf")
ggplot(A3A_B10A_var_log2change_mean, aes(x = wt_mean, y = ko_mean, color = change)) + 
  geom_point(size = 3) + 
  ggtitle("A3A vs B10A") + 
  labs(y = bquote(log[2]~" (KO interactions)"), x = bquote(log[2]~" (WT interactions)")) + 
  scale_x_continuous(expand = c(0, 0), limits = c(1, 4.5)) + 
  scale_y_continuous(expand = c(0, 0), limits = c(1, 4.5)) + 
  geom_abline(intercept = 0, slope = 1, size = 1)
ggplot(A3A_D2A_var_log2change_mean, aes(x = wt_mean, y = ko_mean, color = change)) + 
  geom_point(size = 3) + 
  ggtitle("A3A vs D2A") + 
  labs(y = bquote(log[2]~" (KO interactions)"), x = bquote(log[2]~" (WT interactions)")) + 
  scale_x_continuous(expand = c(0, 0), limits = c(1, 4.5)) + 
  scale_y_continuous(expand = c(0, 0), limits = c(1, 4.5)) + 
  geom_abline(intercept = 0, slope = 1, size = 1)
ggplot(A3A_G91A_var_log2change_mean, aes(x = wt_mean, y = ko_mean, color = change)) + 
  geom_point(size = 3) + 
  ggtitle("A3A vs G91A") + 
  labs(y = bquote(log[2]~" (KO interactions)"), x = bquote(log[2]~" (WT interactions)")) + 
  scale_x_continuous(expand = c(0, 0), limits = c(1, 4.5)) + 
  scale_y_continuous(expand = c(0, 0), limits = c(1, 4.5)) + 
  geom_abline(intercept = 0, slope = 1, size = 1)
ggplot(A3B_B10B_var_log2change_mean, aes(x = wt_mean, y = ko_mean, color = change)) + 
  geom_point(size = 3) + 
  ggtitle("A3B vs B10B") + 
  labs(y = bquote(log[2]~" (KO interactions)"), x = bquote(log[2]~" (WT interactions)")) + 
  scale_x_continuous(expand = c(0, 0), limits = c(1, 4.5)) + 
  scale_y_continuous(expand = c(0, 0), limits = c(1, 4.5)) + 
  geom_abline(intercept = 0, slope = 1, size = 1)
ggplot(A3B_D2C_var_log2change_mean, aes(x = wt_mean, y = ko_mean, color = change)) + 
  geom_point(size = 3) + 
  ggtitle("A3B vs D2C") + 
  labs(y = bquote(log[2]~" (KO interactions)"), x = bquote(log[2]~" (WT interactions)")) + 
  scale_x_continuous(expand = c(0, 0), limits = c(1, 4.5)) + 
  scale_y_continuous(expand = c(0, 0), limits = c(1, 4.5)) + 
  geom_abline(intercept = 0, slope = 1, size = 1)
ggplot(A3B_G91B_var_log2change_mean, aes(x = wt_mean, y = ko_mean, color = change)) + 
  geom_point(size = 3) + 
  ggtitle("A3B vs G91B") + 
  labs(y = bquote(log[2]~" (KO interactions)"), x = bquote(log[2]~" (WT interactions)")) + 
  scale_x_continuous(expand = c(0, 0), limits = c(1, 4.5)) + 
  scale_y_continuous(expand = c(0, 0), limits = c(1, 4.5)) + 
  geom_abline(intercept = 0, slope = 1, size = 1)
dev.off()

##################################################################################################################################
################################################ DIFFERENTIAL INTERACTION HEATMAP ################################################
##################################################################################################################################

diff <- c("B10", "D2", "G91")
for (i in diff) {
  df <- read.csv(paste0("D:/Deitsch_collab/hicdcplus/diff_analysis/diff_res",i,"overA3.txt"), header = TRUE, sep = "\t")
  assign(paste0(i, "overA3_df"), df)
}

sizes <- read.csv('D:/organism_genome/Pfalciparum3D7/PlasmoDB-50_Pfalciparum3D7_Genome.chrom.sizes', header = FALSE, sep = "\t")[1:14,]
colnames(sizes) <- c("chr", "length")
starts <- data.frame(matrix(NA, nrow = nrow(sizes), ncol = ncol(sizes)))
colnames(starts) <- c("chr", "start_bin")
for (i in 1:nrow(sizes)) {
  starts[i, 1] <- sizes[i, 1]
  starts[1, 2] <- 0
  if (i > nrow(sizes) - 1) {
    break
  } else {
    starts[i+1, 2] <- ceiling(sizes[i, 2] / 10000) + starts[i, 2]
  }
}
sizes$length <- ceiling(sizes$length / 10000)

var <- read.csv('D:/organism_genome/Pfalciparum3D7/var_genes.gff', header = FALSE, sep = "\t")[,c(1,4,5,9)]
colnames(var) <- c("chr", "start", "end", "gene")
for(i in 1:nrow(var)) {
  var[i, 2] <- ceiling(var[i, 2] / 10000)
  var[i, 3] <- ceiling(var[i, 3] / 10000)
}
var$gene <- str_extract_all(var$gene, "(?<==).*(?=;)")
colnames(var)[colnames(var) == "start"] <- "startI"

min_max <- list()
j <- 1
for (i in paste0(diff, "overA3_df")) {
  df = get(i)
  df$startI <- df$startI / 10000 + 1
  df$startJ <- df$startJ / 10000 + 1
  df$startIplus1 <- df$startI + 1
  df$startIplus2 <- df$startI + 2
  df$startJplus2 <- df$startJ + 2
  df$startJplus1 <- df$startJ + 1
  df <- df[!(df$startI == df$startJ),]
  df <- df[!(df$startIplus1 == df$startJ),]
  df <- df[!(df$startIplus2 == df$startJ),]
  df <- df[!(df$startI == df$startJplus1),]
  df <- df[!(df$startI == df$startJplus2),]
  min_max[[j]] <- max(abs(df$log2FoldChange))
  assign(i, df)
  j = j + 1
}

quantile_breaks <- function(xs, n = 20) {
  breaks <- quantile(xs, probs = seq(0, 1, length.out = n))
  breaks[!duplicated(breaks)]
}

pdf("D:/Deitsch_collab/figures/diff_heatmaps/diff_intra_heatmaps.pdf")

for (i in paste0(diff, "overA3_df")) {
  mat_breaks <- c()
  df <- get(i)
  list2env(split(df, df$chr), envir = .GlobalEnv)
  chr_lst <- list(unique(df$chr))
  
  for (j in chr_lst) {
    for (k in j) {
      chr = k
      chr_df = get(k)
      chr_mat <- matrix(0, sizes$length[which(sizes$chr == k)], sizes$length[which(sizes$chr == k)])
      chr_mat[as.matrix(chr_df[c("startI", "startJ")])] <- chr_df[["log2FoldChange"]]
      chr_lmat <- t(chr_mat)[,-ncol(chr_mat)]
      chr_mat[lower.tri(chr_mat)] <- chr_lmat[lower.tri(chr_lmat, diag = FALSE)]

      chr_mat <- apply(t(chr_mat),2,rev)
      rownames(chr_mat) <- 1:nrow(chr_mat)
      colnames(chr_mat) <- 1:ncol(chr_mat)
      
      j <- 1
      for (g in var[which(var$chr == chr),]$startI) {
        rownames(chr_mat)[g] <- var[which(var$chr == chr),]$gene[j]
        colnames(chr_mat)[g] <- var[which(var$chr == chr),]$gene[j]
        j = j + 1
      }
      
      rownames(chr_mat) <- rev(rownames(chr_mat))
      
      mat_breaks <- quantile_breaks(chr_df$log2FoldChange, n = 100)
      
      if (length(mat_breaks[which(mat_breaks < 0)]) > length(mat_breaks[which(mat_breaks > 0)])) {
        mat_breaks <- mat_breaks[which(mat_breaks < 0)]
        mat_breaks <- c(-max(abs(chr_df$log2FoldChange)), mat_breaks, 0, -rev(mat_breaks), max(abs(chr_df$log2FoldChange)))
      } else {
        mat_breaks <- mat_breaks[which(mat_breaks > 0)]
        mat_breaks <- c(-max(abs(chr_df$log2FoldChange)), -rev(mat_breaks), 0, mat_breaks, max(abs(chr_df$log2FoldChange)))
      }
      
      if (max(abs(chr_df$log2FoldChange)) == max(abs(mat_breaks))) {
        mat_breaks <- head(mat_breaks, -1)
        mat_breaks <- tail(mat_breaks, -1)
      }
      mat_color <- colorRampPalette(c("blue", "white", "red"))(length(mat_breaks))
      
      annoRow <- data.frame(varGenes = factor(rownames(chr_mat) %in% var$gene))
      rownames(annoRow) <- rownames(chr_mat)
      ann_colors = list(varGenes = c("FALSE"="white", "TRUE"="green"))
      
      pheatmap(chr_mat, 
               main = paste0(gsub("_df", "", i), " ", k), 
               show_colnames = FALSE, 
               show_rownames = FALSE, 
               cluster_cols = FALSE, 
               cluster_rows = FALSE, 
               annotation_row = annoRow, 
               annotation_col = annoRow, 
               annotation_colors = ann_colors, 
               color = mat_color, 
               breaks = mat_breaks, 
               cellwidth = 360 / ncol(chr_mat), 
               cellheight = 360 / ncol(chr_mat), 
               border_color = FALSE)
      
      mat_breaks <- c()
    }
  }
}
dev.off()

######################################### IF YOU WANT THE QUANTILE BREAKS CENTERED ON 0 ##########################################
mat_breaks <- quantile_breaks(chr_df$log2FoldChange, n = 20)

if (length(mat_breaks[which(mat_breaks < 0)]) > length(mat_breaks[which(mat_breaks > 0)])) {
  mat_breaks <- mat_breaks[which(mat_breaks < 0)]
} else {
  mat_breaks <- mat_breaks[which(mat_breaks > 0)]
}

mat_breaks <- c(-max(abs(chr_df$log2FoldChange)), -rev(mat_breaks), 0, mat_breaks, max(abs(chr_df$log2FoldChange)))

if (max(abs(chr_df$log2FoldChange)) == max(abs(mat_breaks))) {
  mat_breaks <- head(mat_breaks, -1)
  mat_breaks <- tail(mat_breaks, -1)
}
##################################################################################################################################

####################################### HiCDC+ doesn't find interchromosomal interactions? #######################################

pdf("D:/Deitsch_collab/figures/diff_heatmaps/diff_inter_heatmap.pdf")

for (i in paste0(diff, "overA3_df")) {
  df = get(i)
  df$startI <- df$startI / 10000 + 1
  df$startJ <- df$startJ / 10000 + 1
  
  mat <- matrix(0, max(max(df$startI), max(df$startJ)), max(max(df$startI), max(df$startJ)))
  mat[as.matrix(df[c("startI", "startJ")])] <- df[["log2FoldChange"]]
  lmat <- t(mat)[,-ncol(mat)]
  mat[lower.tri(mat)] <- lmat[lower.tri(lmat, diag = FALSE)]
  diag(mat) <- 0
  row.names(mat) <- 1:nrow(mat)
  colnames(mat) <- 1:ncol(mat)
  
  breaks <- c(seq(-unlist(max_min[which.max(max_min)]), 0, length.out = ceiling(paletteLength / 2) + 1), 
                  seq(unlist(max_min[which.max(max_min)]) / paletteLength, unlist(max_min[which.max(max_min)]), 
                      length.out = floor(paletteLength / 2)))
  
  for (l in as.list(df$log2FoldChange)) {
    ind <- which(df$log2FoldChange == l, arr.ind = TRUE)
    sI <- df$startI[ind]
    gene <- unlist(df$chr[ind])
    if (!is.null(gene)) {
      colnames(mat)[sI] <- gene[1]
      rownames(mat)[sI] <- gene[1]
    }
  }
  
  pheatmap(mat, 
           main = i, 
           show_colnames = FALSE, 
           show_rownames = FALSE, 
           cluster_cols = FALSE, 
           cluster_rows = FALSE, 
           annotation_row = annoRow, 
           annotation_col = annoRow, 
           #annotation_colors = annoCol, 
           color = mat_color, 
           breaks = breaks, 
           cellwidth = 350 / ncol(mat), 
           cellheight = 350 / ncol(mat), 
           border_color = NA)
}
dev.off()