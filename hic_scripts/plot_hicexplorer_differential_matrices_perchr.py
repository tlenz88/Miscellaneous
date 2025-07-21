#!/usr/bin/env python3

import sys
import argparse
import os
import subprocess
import tempfile
import pandas as pd
from reportlab.pdfgen import canvas
from reportlab.lib.pagesizes import letter, landscape
from PyPDF2 import PdfReader, PdfWriter, PdfFileMerger, PdfFileWriter
import math
import h5py
import numpy as np
import scipy.sparse
import hdf5plugin

def read_chrom_sizes(chrom_sizes_file):
    """Read chromosome sizes from a .chrom.sizes file."""
    sizes = {}
    with open(chrom_sizes_file, 'r') as f:
        for line in f:
            chrom, size = line.strip().split()
            sizes[chrom] = int(size)
    return sizes

def find_data_range(h5_file, chrom_start_bin, start, end, resolution):
    with h5py.File(h5_file, "r") as h5_file:
        # Access the matrix components
        try:
            data = h5_file["matrix/data"][:]
            indices = h5_file["matrix/indices"][:]
            indptr = h5_file["matrix/indptr"][:]
            shape = h5_file["matrix/shape"][:]
        except KeyError:
            print("Error: Matrix components not found. Check the dataset paths.")
            sys.exit()

        # Reconstruct the sparse matrix using the components
        sparse_matrix = scipy.sparse.csr_matrix((data, indices, indptr), shape=shape)
        
        # Find the bin indices corresponding to start and end coordinates
        start_bin = start // resolution #+ chrom_start_bin
        end_bin = end // resolution #+ chrom_start_bin
        
        # Extract the region of interest
        region_values = sparse_matrix[start_bin:end_bin, start_bin:end_bin].toarray()

    # Calculate the min and max values for this region
    vmin = np.min(region_values)
    vmax = np.max(region_values)
    return vmin, vmax


def create_region_plot(h5_file, region, vmin, vmax, output_file, font_size, compartments, tads, loops):
    """Create a plot for a specific region using hicPlotMatrix."""
    cmd = [
        'hicPlotMatrix',
        '--matrix', h5_file,
        '--region', f"{region['chrom']}:{region['start']}-{region['end']}",
        '--outFileName', output_file,
        #'--vMin', str(vmin),
        #'--vMax', str(vmax),
        '--tads', tads,
        '--bigwig', compartments,
        '--loops', loops,
        '--colorMap', 'afmhot_r',
        '--log1p'
    ]
    subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

def create_diff_region_plot(h5_file, region, vmin, vmax, output_file, font_size):
    """Create a plot for a specific region using hicPlotMatrix."""
    cmd = [
        'hicPlotMatrix',
        '--matrix', h5_file,
        '--region', f"{region['chrom']}:{region['start']}-{region['end']}",
        '--outFileName', output_file,
        '--vMin', str(vmin),
        '--vMax', str(vmax),
        '--colorMap', 'RdBu_r'
    ]
    subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

def create_multi_plot_page(plot_files, output_file):
    """Create a single page with multiple plots arranged horizontally."""
    # Read all input PDFs
    readers = [PdfReader(plot_file) for plot_file in plot_files]

    # Get dimensions of first plot (assuming all plots have same dimensions)
    plot_width = float(readers[0].pages[0].mediabox[2])
    plot_height = float(readers[0].pages[0].mediabox[3])

    # Calculate dimensions for the output page
    # Use landscape orientation
    page_width = plot_width * len(plot_files)
    page_height = plot_height

    # Create a new PDF with calculated dimensions
    output = PdfWriter()
    page = output.add_blank_page(width=page_width, height=page_height)

    # Add each plot to the page
    for i, reader in enumerate(readers):
        plot = reader.pages[0]
        plot.scale_to(width=plot_width, height=plot_height)
        page.mergeTranslatedPage(plot, tx=i*plot_width, ty=0)

    # Save the merged page
    with open(output_file, 'wb') as f:
        output.write(f)

def main():
    parser = argparse.ArgumentParser(description='Generate HiC plots for multiple regions')
    #parser.add_argument('--h5_files', nargs='+', required=True, help='Input .h5 files')
    #parser.add_argument('--diff_h5_file', required=True, help='Differential .h5 file')
    parser.add_argument('--samples', nargs='+', required=True, help='Sample names')
    parser.add_argument('--window_size', type=int, required=True, help='Window size for scanning matrices')
    parser.add_argument('--resolution', type=int, required=True, default=50000)
    parser.add_argument('--chrom_sizes', required=True, help='Chromosome sizes file')
    parser.add_argument('--font_size', type=int, default=12, help='Font size for plots')
    parser.add_argument('--bigwig', nargs='+', help='Bigwig files')
    args = parser.parse_args()

    sample_names = args.samples
    resolution = args.resolution

    # Create temporary directory for individual plots
    with tempfile.TemporaryDirectory() as temp_dir:
        # Process whole chromosomes
        chrom_sizes = read_chrom_sizes(args.chrom_sizes)
        chrom_start_bin = 0

        for chrom, size in chrom_sizes.items():
            final_pages = []
            chrom_start_bin += size // resolution + 1

            end_counter = 0
            for start in list(range(0, int(size), int(args.window_size)//2)):
                end = start + args.window_size
                if end > size:
                    end = size
                    end_counter += 1
                    if end_counter > 1:
                        break
                
                region = {
                    'chrom': chrom,
                    'start': start,
                    'end': end,
                    'name': chrom
                }

                max_values = []
                for idx, sample_name in enumerate(sample_names):
                    input_file = f'/mnt/f/toxo_project/HiC/Hsapien_output/output_files/{sample_name}_HiC/hicexplorer_files_two_reps/{resolution}/{sample_name}_{resolution}_{chrom}_corrected.h5'
                    max_values.extend(find_data_range(input_file, chrom_start_bin, start, end, resolution))
                vmax = max(max_values)
                
                # Create plots for each input file
                print(f'Plotting {chrom}:{start}-{end}')
                temp_pdfs = []
                for idx, sample_name in enumerate(sample_names):
                    output_file = os.path.join(temp_dir, f"{chrom}_{start}_{idx}.pdf")
                    input_file = f'/mnt/f/toxo_project/HiC/Hsapien_output/output_files/{sample_name}_HiC/hicexplorer_files_two_reps/{resolution}/{sample_name}_{resolution}_{chrom}_corrected.h5'
                    compartments = f'/mnt/f/toxo_project/HiC/Hsapien_output/output_files/{sample_name}_HiC/hicexplorer_files_two_reps/{resolution}/{sample_name}_{resolution}_{chrom}_PCA1.bw'
                    tads = f'/mnt/f/toxo_project/HiC/Hsapien_output/output_files/{sample_name}_HiC/hicexplorer_files_two_reps/{resolution}/{sample_name}_{resolution}_{chrom}_domains.bed'
                    loops = f'/mnt/f/toxo_project/HiC/Hsapien_output/output_files/{sample_name}_HiC/hicexplorer_files_two_reps/{resolution}/{sample_name}_{resolution}_{chrom}_loops.bedgraph'
                    create_region_plot(input_file, region, 0, vmax, output_file, args.font_size, compartments, tads, loops)
                    temp_pdfs.append(output_file)
                
                for sample1, sample2 in zip(sample_names[0::2], sample_names[1::2]):
                    diff_h5_file = f'/mnt/f/toxo_project/HiC/Hsapien_output/output_files/{sample2}_HiC/hicexplorer_files_two_reps/{resolution}/{sample2}_vs_{sample1}_{resolution}_{chrom}_normalized.h5'
                    vmin, vmax = find_data_range(diff_h5_file, chrom_start_bin, start, end, resolution)
                    vmax = max(-vmin, vmax)

                    output_file = os.path.join(temp_dir, f"{chrom}_{start}_diff.pdf")
                    create_diff_region_plot(diff_h5_file, region, -vmax, vmax, output_file, args.font_size)
                    temp_pdfs.append(output_file)
                
                # Create single page with all plots
                page_output = os.path.join(temp_dir, f"{chrom}_{start}_merged.pdf")
                create_multi_plot_page(temp_pdfs, page_output)
                final_pages.append(page_output)

            # Combine all pages into final PDF
            merger = PdfFileMerger()
            for page_file in final_pages:
                merger.append(page_file)

            for sample1, sample2 in zip(sample_names[0::2], sample_names[1::2]):
                window_str = int(args.window_size / 1e6)
                outfile = f'/mnt/f/toxo_project/HiC/Hsapien_output/output_files/{sample2}_HiC/hicexplorer_files_two_reps/{resolution}/{sample2}_vs_{sample1}_{resolution}_{chrom}_{window_str}Mb.pdf'
                # Save final PDF
                with open(outfile, 'wb') as f:
                    merger.write(f)
                    merger.close()

if __name__ == "__main__":
    main()
