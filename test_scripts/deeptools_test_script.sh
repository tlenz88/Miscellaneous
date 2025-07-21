#!/usr/bin/env bash

# Merge replicates for experimental and input samples of condition 1
samtools merge -@ 18 -@ 4 -o condition1_H3K9_merged.bam condition1_H3K9_rep1.bam condition1_H3K9_rep2.bam
samtools merge -@ 18 -@ 4 -o condition1_input_merged.bam condition1_input_rep1.bam condition1_input_rep2.bam

# Merge replicates for experimental and input samples of condition 2
samtools merge -@ 18 -@ 4 -o condition2_H3K9_merged.bam condition2_H3K9_rep1.bam condition2_H3K9_rep2.bam
samtools merge -@ 18 -@ 4 -o condition2_input_merged.bam condition2_input_rep1.bam condition2_input_rep2.bam

# Create index for condition1
samtools index -@ 18 condition1_H3K9_merged.bam
samtools index -@ 18 condition1_input_merged.bam

# Create index for condition2
samtools index -@ 18 condition2_H3K9_merged.bam
samtools index -@ 18 condition2_input_merged.bam

# Subtract input from experimental for each condition
bamCompare -b1 condition1_H3K9_merged.bam -b2 condition1_input_merged.bam --outFileName condition1_input_sub.bw --scaleFactorsMethod None --normalizeUsing RPKM --operation subtract -bs 1 --ignoreDuplicates -p 4
bamCompare -b1 condition2_H3K9_merged.bam -b2 condition2_input_merged.bam --outFileName condition2_input_sub.bw --scaleFactorsMethod None --normalizeUsing RPKM --operation subtract -bs 1 --ignoreDuplicates -p 4

# Find TSS coverage for both conditions
computeMatrix reference-point -S condition1_input_sub.bw condition2_input_sub.bw -R genes.bed -a 500 -b 1500 -o TSS_coverage.matrix.gz -p 4 -bs 1

# Create profile plot of TSS coverage
plotProfile -m TSS_coverage.matrix.gz -o TSS_coverage_profile.png --perGroup --colors red blue --plotTitle "H3K9 Enrichment" --samplesLabel "condition1" "condition2"
