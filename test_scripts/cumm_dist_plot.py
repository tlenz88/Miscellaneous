import pybedtools
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats

def read_bed_file(file_path):
    return pybedtools.BedTool(file_path)

def cumulative_coverage(peaks, genes):
    coverages = []
    for gene in genes:
        overlap = sum(min(peak.end, gene.end) - max(peak.start, gene.start) 
                      for peak in peaks if peak.chrom == gene.chrom and peak.start < gene.end and peak.end > gene.start)
        coverages.append(overlap / (gene.end - gene.start))
    return sorted(coverages)

# Read the peak files
ko_peaks = read_bed_file('D2_H3K9me3_peaks.broadPeak')
wt_peaks = read_bed_file('A3_H3K9me3_peaks.broadPeak')

# Read the var and rifin genes file
var_rifin_genes = read_bed_file('/mnt/f/organism_genome/Pfalciparum3D7/var.bed')

# Calculate cumulative coverage
ko_coverage = cumulative_coverage(ko_peaks, var_rifin_genes)
wt_coverage = cumulative_coverage(wt_peaks, var_rifin_genes)

# Create the plot
plt.figure(figsize=(10, 6))
plt.plot(ko_coverage, np.linspace(0, 1, len(ko_coverage)), label='KO')
plt.plot(wt_coverage, np.linspace(0, 1, len(wt_coverage)), label='WT')
plt.xlabel('Fraction of gene covered by H3K9me3')
plt.ylabel('Cumulative fraction of genes')
plt.legend()
plt.title('Cumulative H3K9me3 Coverage of var and rifin Genes')
plt.savefig('cumulative_coverage_plot.png')
plt.close()

# Perform statistical test
_, p_value = stats.ks_2samp(ko_coverage, wt_coverage)
print(f"Kolmogorov-Smirnov test p-value: {p_value}")
