#!/usr/bin/env python3

import pandas as pd
import numpy as np
import sys
import os
import subprocess
from pathlib import Path
import tempfile
from PyPDF2 import PdfReader, PdfWriter, PdfFileMerger, PdfFileWriter
import configparser

def read_bed_file(file_path, file_type):
    """
    Read BED files with specific column names based on file type
    """
    if file_type == 'tad':
        columns = ['chromosome', 'start', 'end', 'name', 'score', 'strand', 
                  'pval_left', 'pval_right', 'pval_intra', 
                  'w_left', 'w_right', 'w_intra']
        # Skip commented header lines starting with #
        df = pd.read_csv(file_path, sep='\t', names=columns, comment='#')
    elif file_type == 'deg':
        columns = ['chromosome', 'start', 'end', 'gene', 'log2fc', 'strand']
        df = pd.read_csv(file_path, sep='\t', names=columns)
    elif file_type == 'domains':
        columns = ['chromosome', 'start', 'end', 'name', 'score', 'strand', 
                  'start2', 'end2', 'other']
        df = pd.read_csv(file_path, sep='\t', names=columns)
    
    return df

def is_overlapping(row1, row2):
    if (row1['start'] < row2['end']) and (row1['end'] > row2['start']):
        return True
    
    return False

def filter_overlapping_rows(df1, df2):
    # List to store the rows from df1 that overlap with df2
    overlapping_rows = []

    # Iterate over each row in df1
    for idx1, row1 in df1.iterrows():
        overlap_found = False
        # Iterate over each row in df2 to check if there's any overlap
        for idx2, row2 in df2.iterrows():
            if is_overlapping(row1, row2):
                overlap_found = True
                break  # If overlap found, no need to check further
        
        # If overlap is found, keep the row from df1
        if overlap_found:
            overlapping_rows.append(row1)

    # Return a new dataframe with only overlapping rows
    return pd.DataFrame(overlapping_rows)

def find_containing_tad(gene_pos, domains_df):
    """
    Find the TAD that contains a given gene position and return its index
    """
    chr_domains = domains_df[domains_df['chromosome'] == gene_pos['chromosome']].sort_values('start')
    containing_tad = chr_domains[
        (chr_domains['start'] <= gene_pos['start']) &
        (chr_domains['end'] >= gene_pos['end'])
    ]
    if containing_tad.empty:
        return None, None
    return containing_tad.iloc[0], chr_domains.index.get_loc(containing_tad.index[0])

def find_gene_cluster(gene_idx, genes_df, domains_df):
    """
    Find cluster of genes that are in adjacent TADs starting from the given gene
    Returns list of (gene, tad, tad_index) tuples and the last downstream TAD index
    """
    chromosome = genes_df.iloc[gene_idx]['chromosome']
    chr_domains = domains_df[domains_df['chromosome'] == chromosome].sort_values('start')
    
    cluster = []
    current_gene_idx = gene_idx
    last_tad_idx = None
    
    while current_gene_idx < len(genes_df):
        current_gene = genes_df.iloc[current_gene_idx]
        
        # Only process genes from the same chromosome
        if current_gene['chromosome'] != chromosome:
            break
            
        # Find TAD containing this gene
        tad, tad_idx = find_containing_tad(current_gene, domains_df)
        
        if tad is None:
            current_gene_idx += 1
            continue
            
        # If this is the first gene or the TAD is adjacent to the last one
        if not cluster or tad_idx <= last_tad_idx + 1:
            cluster.append((current_gene, tad, tad_idx))
            last_tad_idx = tad_idx
            current_gene_idx += 1
        else:
            break
    
    # Get the downstream TAD index (one past the last gene's TAD)
    downstream_tad_idx = last_tad_idx + 1 if last_tad_idx is not None else None
    
    return cluster, downstream_tad_idx

def get_expanded_region(gene_cluster, downstream_tad_idx, domains_df):
    """
    Get the expanded region coordinates including upstream TAD of first gene
    and downstream TAD of last gene
    """
    if not gene_cluster:
        return None
    
    chromosome = gene_cluster[0][0]['chromosome']
    chr_domains = domains_df[domains_df['chromosome'] == chromosome].sort_values('start')
    
    # Get first gene's TAD index and upstream TAD
    first_gene_tad_idx = gene_cluster[0][2]
    upstream_tad = chr_domains.iloc[first_gene_tad_idx - 1] if first_gene_tad_idx > 0 else None
    
    # Get downstream TAD after the last gene's TAD
    downstream_tad = chr_domains.iloc[downstream_tad_idx] if downstream_tad_idx < len(chr_domains) else None
    
    # Define region boundaries
    start = upstream_tad['start'] if upstream_tad is not None else gene_cluster[0][1]['start']
    end = downstream_tad['end'] if downstream_tad is not None else gene_cluster[-1][1]['end']
    
    return {
        'chromosome': chromosome,
        'start': start,
        'end': end
    }

def create_plot_command(region, ini_file, output_file):
    """
    Create hicPlotTADs command for the given region
    """
    cmd = [
        'hicPlotTADs',
        '--tracks', ini_file,
        '--region', f"{region['chromosome']}:{region['start']}-{region['end']}",
        '--outFileName', output_file,
    ]
    return cmd

def main():
    if len(sys.argv) != 8:
        print("Usage: script.py <deg_file> <diff_tads_file> <domains.bed> <tracks.ini> <chromosome>")
        sys.exit(1)

    # Read input files
    deg_file = sys.argv[1]
    diff_tads_file = sys.argv[2]
    domains_file = sys.argv[3]
    ini_file = sys.argv[4]
    chromosome = sys.argv[5]
    resolution = sys.argv[6]
    condition = sys.argv[7]

    # Create output directory
    output_dir = Path('tad_plots')
    output_dir.mkdir(exist_ok=True)

    # Read the files
    try:
        deg_df = read_bed_file(deg_file, 'deg')
        print(f"Successfully loaded {len(deg_df)} differentially expressed genes")
        
        diff_tads_df = read_bed_file(diff_tads_file, 'tad')
        print(f"Successfully loaded {len(diff_tads_df)} differential TADs")
        
        domains_df = read_bed_file(domains_file, 'domains')
        print(f"Successfully loaded {len(domains_df)} domains")
        
    except Exception as e:
        print(f"Error reading input files: {e}")
        sys.exit(1)


    # Sort DEGs by chromosome and start position
    deg_df = deg_df.sort_values(['chromosome', 'start']).reset_index(drop=True)
    deg_df = deg_df[deg_df['chromosome'] == chromosome]
    deg_df = filter_overlapping_rows(deg_df, diff_tads_df)

    # Process genes
    gene_idx = 0
    plots_generated = 0
    
    with tempfile.TemporaryDirectory() as temp_dir:
        final_pages = []

        while gene_idx < len(deg_df):
            try:
                # Find cluster of genes in adjacent TADs
                gene_cluster, downstream_tad_idx = find_gene_cluster(gene_idx, deg_df, domains_df)

                if gene_cluster:
                    # Get expanded region for the cluster
                    region = get_expanded_region(gene_cluster, downstream_tad_idx, domains_df)

                    if region:
                        # Create descriptive filename using first and last gene in cluster
                        plot_coords = f"{region['chromosome']}_{region['start']}_to_{region['end']}"
                        first_gene = gene_cluster[0][0]['gene']
                        last_gene = gene_cluster[-1][0]['gene']
                        if first_gene == last_gene:
                            filename = f"{condition}_vs_HFF_{resolution}_{plot_coords}_{first_gene}_diff_TADs.pdf"
                        else:
                            filename = f"{condition}_vs_HFF_{resolution}_{plot_coords}_{first_gene}_to_{last_gene}_diff_TADs.pdf"
                        
                        output_file = output_dir / filename
                        final_pages.append(output_file)

                        # Generate plot
                        plot_cmd = create_plot_command(region, ini_file, str(output_file))
                        subprocess.run(plot_cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
                        plots_generated += 1
                        
                        # Print information about the cluster
                        print(f"\nGenerated plot for gene cluster at {output_file}")
                        print(f"Genes in cluster:")
                        for gene, tad, _ in gene_cluster:
                            print(f"  Gene: {gene['gene']}, Log2FC: {gene['log2fc']}")
                            print(f"  TAD position: {tad['chromosome']}:{tad['start']}-{tad['end']}")
                        print(f"Plot region: {region['chromosome']}:{region['start']}-{region['end']}")
                        print("-" * 50)
                    
                    # Move to the next gene after the cluster
                    gene_idx += len(gene_cluster)
                else:
                    gene_idx += 1
                    
            except Exception as e:
                print(f"Error processing gene at index {gene_idx}: {e}")
                gene_idx += 1
                continue

        # Combine all pages into final PDF
        merger = PdfFileMerger()
        for page_file in final_pages:
            merger.append(page_file)

        outfile = os.path.join(output_dir, f"{condition}_vs_HFF_{resolution}_{chromosome}_TAD_DEGenes.pdf")

        # Save final PDF
        with open(outfile, 'wb') as f:
            merger.write(f)
            merger.close()

    # Print summary statistics
    print("\nSummary:")
    print(f"Total differential genes processed: {len(deg_df)}")
    print(f"Successful plots generated: {plots_generated}")

if __name__ == "__main__":
    main()
