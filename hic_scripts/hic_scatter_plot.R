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
for(i in 1:nrow(var)) {
  var[i, 5] <- ceiling(((var[i, 3] - var[i, 2]) / 2 + var[i, 2]) / 10000)
}
colnames(var) <- c("chr", "start", "end", "gene", "bin")
var$gene <- str_extract_all(var$gene, "(?<==).*(?=;)")
var <- var %>% mutate_at("gene", as.character)
var_bins <- var$bin

sizes <- read.csv('F:/organism_genome/Pfalciparum3D7/PlasmoDB-50_Pfalciparum3D7_Genome.chrom.sizes', header = FALSE, sep = "\t")[1:14,]
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
var["bin"] <- var["bin"] + var["start_bin"]
var <- subset(var, select = -start_bin)
var <- subset(var, select = c("chr", "gene", "bin"))
var <- var[!is.na(var$bin),]
var_bins <- var$bin

koA <- c("AP2_16A_KO", "AP2_40A_KO")
for (i in koA) {
  df = read.csv(paste0("F:/ap2_hic/output_files/",i,"_10000_iced_cpm.matrix"), header = FALSE, sep = "\t")
  colnames(df) = c("bin1", "bin2", "ko_intA")
  assign(paste0(i), df)
}

koB <- c("AP2_16B_KO", "AP2_40B_KO")
for (i in koB) {
  df = read.csv(paste0("F:/ap2_hic/output_files/",i,"_10000_iced_cpm.matrix"), header = FALSE, sep = "\t")
  colnames(df) = c("bin1", "bin2", "ko_intB")
  assign(paste0(i), df)
}

AP2_16_KO <- left_join(AP2_16A_KO, AP2_16B_KO, by=c("bin1", "bin2"))
AP2_16_KO[is.na(AP2_16_KO)] = 0
AP2_16_KO$ko_int <- AP2_16_KO$ko_intA * (sum(AP2_16_KO$ko_intA) / (sum(AP2_16_KO$ko_intA) + sum(AP2_16_KO$ko_intB))) + 
  (AP2_16_KO$ko_intB * (sum(AP2_16_KO$ko_intB) / (sum(AP2_16_KO$ko_intA) + sum(AP2_16_KO$ko_intB))))
AP2_16_KO$ko_int <- AP2_16_KO$ko_int / sum(AP2_16_KO$ko_int) * 1000000

AP2_40_KO <- left_join(AP2_40A_KO, AP2_40B_KO, by=c("bin1", "bin2"))
AP2_40_KO[is.na(AP2_40_KO)] = 0
AP2_40_KO$ko_int <- AP2_40_KO$ko_intA * (sum(AP2_40_KO$ko_intA) / (sum(AP2_40_KO$ko_intA) + sum(AP2_40_KO$ko_intB))) + 
  (AP2_40_KO$ko_intB * (sum(AP2_40_KO$ko_intB) / (sum(AP2_40_KO$ko_intA) + sum(AP2_40_KO$ko_intB))))
AP2_40_KO$ko_int <- AP2_40_KO$ko_int / sum(AP2_40_KO$ko_int) * 1000000

for (i in c("AP2_16A_CTRL", "AP2_40A_CTRL")) {
  df = read.csv(paste0("F:/ap2_hic/output_files/",i,"_10000_iced_cpm.matrix"), header = FALSE, sep = "\t")
  colnames(df) = c("bin1", "bin2", "wt_intA")
  assign(paste0(i), df)
}

for (i in c("AP2_16B_CTRL", "AP2_40B_CTRL")) {
  df = read.csv(paste0("F:/ap2_hic/output_files/",i,"_10000_iced_cpm.matrix"), header = FALSE, sep = "\t")
  colnames(df) = c("bin1", "bin2", "wt_intB")
  assign(paste0(i), df)
}

AP2_16_wt <- left_join(AP2_16A_CTRL, AP2_16B_CTRL, by=c("bin1", "bin2"))
AP2_16_wt[is.na(AP2_16_wt)] = 0
AP2_16_wt$wt_int <- AP2_16_wt$wt_intA * (sum(AP2_16_wt$wt_intA) / (sum(AP2_16_wt$wt_intA) + sum(AP2_16_wt$wt_intB))) + 
  (AP2_16_wt$wt_intB * (sum(AP2_16_wt$wt_intB) / (sum(AP2_16_wt$wt_intA) + sum(AP2_16_wt$wt_intB))))
AP2_16_wt$wt_int <- AP2_16_wt$wt_int / sum(AP2_16_wt$wt_int) * 1000000

AP2_40_wt <- left_join(AP2_40A_CTRL, AP2_40B_CTRL, by=c("bin1", "bin2"))
AP2_40_wt[is.na(AP2_40_wt)] = 0
AP2_40_wt$wt_int <- AP2_40_wt$wt_intA * (sum(AP2_40_wt$wt_intA) / (sum(AP2_40_wt$wt_intA) + sum(AP2_40_wt$wt_intB))) + 
  (AP2_40_wt$wt_intB * (sum(AP2_40_wt$wt_intB) / (sum(AP2_40_wt$wt_intA) + sum(AP2_40_wt$wt_intB))))
AP2_40_wt$wt_int <- AP2_40_wt$wt_int / sum(AP2_40_wt$wt_int) * 1000000

for (i in c("AP2_16_KO", "AP2_40_KO", "AP2_16_wt", "AP2_40_wt")) {
  df <- get(i)
  df <- subset(df, bin1 %in% var_bins)
  df <- subset(df, bin2 %in% var_bins)
  assign(paste0(i, "_var"), df)
}

AP2_16_var <- left_join(AP2_16_wt_var, AP2_16_KO_var, by=c("bin1", "bin2"))
AP2_40_var <- left_join(AP2_40_wt_var, AP2_40_KO_var, by=c("bin1", "bin2"))

rm(AP2_16A_CTRL, AP2_16B_CTRL, AP2_16A_KO, AP2_16B_KO, AP2_40A_CTRL, AP2_40B_CTRL, AP2_40A_KO, AP2_40B_KO, 
   AP2_16_KO_var, AP2_16_KO, AP2_16_wt_var, AP2_16_wt, AP2_40_wt_var, AP2_40_wt, AP2_40_KO_var, AP2_40_KO, 
   i, koA, koB, df, sizes, starts, var, var_bins)

for (i in c("AP2_16_var", "AP2_40_var")) {
  df = get(i)
  df[is.na(df)] = 0
  df$change <- df$ko_int > df$wt_int
  df["change"][df["change"] == "TRUE"] <- "up"
  df["change"][df["change"] == "FALSE"] <- "down"
  df$change <- factor(df$change, levels = c("up", "down"))
  assign(i, df)
}

top2.5percent <- max(min(head(sort(AP2_16_var$wt_int, decreasing=TRUE), n=ceiling(length(AP2_16_var$wt_int)*.025))),
                     min(head(sort(AP2_16_var$ko_int, decreasing=TRUE), n=ceiling(length(AP2_16_var$ko_int)*.025))),
                     min(head(sort(AP2_40_var$wt_int, decreasing=TRUE), n=ceiling(length(AP2_40_var$wt_int)*.025))),
                     min(head(sort(AP2_40_var$ko_int, decreasing=TRUE), n=ceiling(length(AP2_40_var$ko_int)*.025))))

twoSD <- max(2*sd(AP2_16_var$wt_int),
             2*sd(AP2_16_var$ko_int),
             2*sd(AP2_40_var$wt_int),
             2*sd(AP2_40_var$ko_int))

ggplot(AP2_16_var, aes(x = wt_int, y = ko_int, color = change)) + 
  geom_point(size = 2) + 
  labs(y = "KO interactions", x = "WT interactions", caption = paste0("Increased var interactions: ", length(which(AP2_16_var$change=="up")), 
                                                                      "  Decreased var interactions: ", length(which(AP2_16_var$change=="down")))) + 
  scale_x_continuous(expand = c(0, 0), limits = c(-.02*twoSD, twoSD)) + 
  scale_y_continuous(expand = c(0, 0), limits = c(-.02*twoSD, twoSD)) + 
  geom_abline(intercept = 0, slope = 1, size = 1)
