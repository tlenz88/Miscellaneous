import numpy as np
import pandas as pd
from diffhic.main import diffHiC
from diffhic.utils import read_hic_matrix

# Load Hi-C data
condition1_reps = ['/mnt/f/toxo_project/HiC/Hsapien_output/HFF_A_HiC/hicexplorer_files/1000000/HFF_A_HiC_corrected.h5', '/mnt/f/toxo_project/HiC/Hsapien_output/HFF_B_HiC/hicexplorer_files/1000000/HFF_B_HiC_corrected.h5', '/mnt/f/toxo_project/HiC/Hsapien_output/HFF_C_HiC/hicexplorer_files/1000000/HFF_C_HiC_corrected.h5']
condition2_reps = ['/mnt/f/toxo_project/HiC/Hsapien_output/ME49RFP_6hpi_A_HiC/hicexplorer_files/1000000/ME49RFP_6hpi_A_HiC_corrected.h5', '/mnt/f/toxo_project/HiC/Hsapien_output/ME49RFP_6hpi_B_HiC/hicexplorer_files/1000000/ME49RFP_6hpi_B_HiC_corrected.h5', '/mnt/f/toxo_project/HiC/Hsapien_output/ME49RFP_6hpi_C_HiC/hicexplorer_files/1000000/ME49RFP_6hpi_C_HiC_corrected.h5']
condition3_reps = ['/mnt/f/toxo_project/HiC/Hsapien_output/ME49RFP_24hpi_A_HiC/hicexplorer_files/1000000/ME49RFP_24hpi_A_HiC_corrected.h5', '/mnt/f/toxo_project/HiC/Hsapien_output/ME49RFP_24hpi_B_HiC/hicexplorer_files/1000000/ME49RFP_24hpi_B_HiC_corrected.h5', '/mnt/f/toxo_project/HiC/Hsapien_output/ME49RFP_24hpi_C_HiC/hicexplorer_files/1000000/ME49RFP_24hpi_C_HiC_corrected.h5']

# Read Hi-C matrices
condition1_matrices = [read_hic_matrix(file, file_format='hdf5') for file in condition1_reps]
condition2_matrices = [read_hic_matrix(file, file_format='hdf5') for file in condition2_reps]
condition3_matrices = [read_hic_matrix(file, file_format='hdf5') for file in condition3_reps]

# Normalize and correct the Hi-C matrices (optional)
condition1_matrices = [normalize_and_correct(matrix) for matrix in condition1_matrices]
condition2_matrices = [normalize_and_correct(matrix) for matrix in condition2_matrices]
condition3_matrices = [normalize_and_correct(matrix) for matrix in condition3_matrices]

# Perform differential analysis with diffHiC
diff_interactions = diffHiC(
    condition1_matrices, 
    condition2_matrices,
    condition3_matrices,
    fdr_cutoff=0.05,
    min_distance=10000,
    max_distance=1000000
)

# Save results to file
diff_interactions.to_csv('differential_interactions.csv', index=False)
print(f'Differential interactions saved to differential_interactions.csv')