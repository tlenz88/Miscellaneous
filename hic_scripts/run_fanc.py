import fanc
from fanc.pairs import generate_pairs_split as generate_pairs
import genomic_regions as gr
from fanc.hic import Hic
from fanc.hic import LowCoverageFilter
from fanc.hic import kr_balancing
import matplotlib.pyplot as plt
import fanc.plotting as fancplot

input_dir = '/mnt/f/toxo_project/HiC/Hsapien_output/output_files'
res = 250000
genome_fa = '/mnt/f/organism_genome/Hsapien_autosomes/GRCh38.fasta'

HFF_hic = fanc.load(f'{input_dir}/HFF_HiC/hicexplorer_files_two_reps/{res}/HFF_{res}.cool')
ME49RFP_6hpi_hic = fanc.load(f'{input_dir}/ME49RFP_6hpi_HiC/hicexplorer_files_two_reps/{res}/ME49RFP_6hpi_{res}.cool')
ME49RFP_24hpi_hic = fanc.load(f'{input_dir}/ME49RFP_24hpi_HiC/hicexplorer_files_two_reps/{res}/ME49RFP_24hpi_{res}.cool')
"""
#m = hic.matrix()
lc_filter = LowCoverageFilter(hic, rel_cutoff=0.2)
hic.filter(lc_filter)
hic.run_queued_filters()
kr_balancing(hic, whole_matrix=False, restore_coverage=False)
intra_expected, intra_expected_chromosome, inter_expected = hic.expected_values()
"""

HFF_hic_ab = fanc.ABCompartmentMatrix.from_hic(HFF_hic)
ME49RFP_6hpi_hic_ab = fanc.ABCompartmentMatrix.from_hic(ME49RFP_6hpi_hic)
ME49RFP_24hpi_hic_ab = fanc.ABCompartmentMatrix.from_hic(ME49RFP_24hpi_hic)

pdf_name = "saddle_plot.pdf"
fig, (ax1, ax2, ax3) = plt.subplots(1, 3)
HFF_hic_mp = fancplot.SquareMatrixPlot(HFF_hic_ab, ax=ax1,
                                       norm='lin', colormap='RdBu_r',
                                       vmin=-1, vmax=1,
                                       draw_minor_ticks=False)
ME49RFP_6hpi_hic_ab_mp = fancplot.SquareMatrixPlot(ME49RFP_6hpi_hic_ab, ax=ax2,
                                                   norm='lin', colormap='RdBu_r',
                                                   vmin=-1, vmax=1,
                                                   draw_minor_ticks=False)
ME49RFP_24hpi_hic_ab_mp = fancplot.SquareMatrixPlot(ME49RFP_24hpi_hic_ab, ax=ax3,
                                                   norm='lin', colormap='RdBu_r',
                                                   vmin=-1, vmax=1,
                                                   draw_minor_ticks=False)
pdf.savefig()
plt.close()

HFF_hic_ev = HFF_hic_ab.eigenvector(genome=genome_fa, force=True)
ME49RFP_6hpi_hic_ev = ME49RFP_6hpi_hic_ab.eigenvector(genome=genome_fa, force=True)
ME49RFP_24hpi_hic_ev = ME49RFP_24hpi_hic_ab.eigenvector(genome=genome_fa, force=True)

profile, cutoffs = HFF_hic.enrichment_profile(HFF_hic, genome=genome_fa)
profile, cutoffs = ME49RFP_6hpi_hic.enrichment_profile(ME49RFP_6hpi_hic, genome=genome_fa)
profile, cutoffs = ME49RFP_24hpi_hic.enrichment_profile(ME49RFP_24hpi_hic, genome=genome_fa)

fig, axes = fancplot.saddle_plot(profile, cutoffs)
