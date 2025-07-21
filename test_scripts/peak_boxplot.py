import pybedtools
import matplotlib.pyplot as plt
from scipy import stats
from matplotlib.backends.backend_pdf import PdfPages

def get_peak_widths(peak_file):
    peaks = pybedtools.BedTool(peak_file)
    return [peak.length for peak in peaks]

ko_widths = get_peak_widths('D2_H3K9me3_peaks.broadPeak')
wt_widths = get_peak_widths('A3_H3K9me3_peaks.broadPeak')

pdf = PdfPages('peak_width_boxplot.pdf')

plt.boxplot([ko_widths, wt_widths], labels=['KO', 'WT'])
plt.ylabel('Peak Width (bp)')
plt.title('H3K9me3 Peak Widths')
pdf.savefig()
plt.close()
pdf.close()

_, p_value = stats.mannwhitneyu(ko_widths, wt_widths)
print(f"Mann-Whitney U test p-value: {p_value}")
