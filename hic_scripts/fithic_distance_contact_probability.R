library(ggpubr)

sig1 <- read.csv("F:/Noetzel_collab/HiC/output_files_BGF/WT_subsample/FitHiC.spline_pass1.res10000.significances.txt.gz", header=TRUE, sep="\t")
pass1 <- read.csv("F:/Noetzel_collab/HiC/output_files_BGF/WT_subsample/FitHiC.fithic_pass1.res10000.txt", header=TRUE, sep="\t")
sig2 <- read.csv("F:/Noetzel_collab/HiC/output_files_BGF/Xp/FitHiC.spline_pass1.res10000.significances.txt.gz", header=TRUE, sep="\t")
pass2 <- read.csv("F:/Noetzel_collab/HiC/output_files_BGF/Xp/FitHiC.fithic_pass1.res10000.txt", header=TRUE, sep="\t")

pass1$logdist <- log10(pass1$avgGenomicDist)
pass1$logprob <- log10(pass1$contactProbability)
pass1$sample <- "WT"
pass2$logdist <- log10(pass2$avgGenomicDist)
pass2$logprob <- log10(pass2$contactProbability)
pass2$sample <- "Xp"

pass <- merge(pass1[,c(6,7,8)], pass2[,c(6,7,8)], all=TRUE)

pdf(file.path("F:/Noetzel_collab/HiC/output_files_BGF/distance_contact_probability.pdf"))
ggscatter(pass, x = "logdist", y = "logprob", color = "sample", palette = c("#004D40", "#FFC107"),
          xlab = "Genomic distance (log10 bp)", ylab = "Contact probabilitiy (log10)", 
          add = "reg.line", conf.int = TRUE, cor.coef = TRUE, cor.method = "pearson")
dev.off()
