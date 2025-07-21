import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import heapq
from collections import defaultdict
from scipy import sparse

def parse_arguments():
    parser = argparse.ArgumentParser(description="Generate intrachromosomal HiC heatmaps.")
    parser.add_argument("-m", "--matrices", nargs="+", required=True, help="List of HiC-Pro matrix files")
    parser.add_argument("-c", "--chrom_sizes", required=True, help="Tab-delimited file containing chromosome sizes")
    parser.add_argument("-r", "--resolution", type=int, required=True, help="Binning resolution")
    parser.add_argument("-o", "--output", default="output.pdf", help="Output filename (PDF)")
    parser.add_argument("-g", "--gff", help="GFF file")
    parser.add_argument("-e", "--centromeres", help="Tab-delimited file containing centromere coordinates")
    parser.add_argument("-l", "--gene_list", nargs="+", help="List of genes to extract from the GFF file")
    parser.add_argument("-n", "--normalization", choices=["raw", "KR", "ICE"], default="raw", help="Normalization method")
    return parser.parse_args()

def load_chrom_bins(chrom_sizes_file, resolution):
    chrom_sizes = pd.read_csv(chrom_sizes_file, sep="\t", header=None, names=["chrom", "size"])
    chrom_bins = {}
    current_bin = 0
    for _, row in chrom_sizes.iterrows():
        n_bins = row['size'] // resolution + 1
        chrom_bins[row['chrom']] = (current_bin, current_bin + n_bins)
        current_bin += n_bins
    return chrom_bins, chrom_sizes

def load_matrix(matrix_file, chrom, chrom_bins, resolution):
    start_bin, end_bin = chrom_bins[chrom]
    df = pd.read_csv(matrix_file, sep="\t", header=None, names=["bin1", "bin2", "count"])
    
    # Filter interactions for the specific chromosome
    df = df[(df["bin1"] >= start_bin) & (df["bin1"] < end_bin) & 
            (df["bin2"] >= start_bin) & (df["bin2"] < end_bin)]
    
    # Adjust bin numbers to start from 0 for each chromosome
    df["bin1"] -= start_bin
    df["bin2"] -= start_bin
    
    # Create a full matrix with zeros for missing interactions
    n_bins = end_bin - start_bin
    full_matrix = sparse.csr_matrix((n_bins, n_bins), dtype=float)
    
    # Fill in the observed interactions
    full_matrix[df["bin1"], df["bin2"]] = df["count"]
    full_matrix[df["bin2"], df["bin1"]] = df["count"]  # Ensure symmetry
    
    # Convert back to DataFrame
    full_df = pd.DataFrame({"bin1": full_matrix.nonzero()[0],
                            "bin2": full_matrix.nonzero()[1],
                            "count": full_matrix.data})
    
    return full_df

def kr_normalize(matrix):
    n = matrix.shape[0]
    tol = 1e-6
    delta = 0.1
    x = np.ones(n)

    for _ in range(50):
        z = matrix.dot(x)
        w = 1.0 / (z + delta)
        y = matrix.T.dot(w)
        y = y / np.linalg.norm(y)
        if np.allclose(x, y, rtol=tol):
            break
        x = y

    return sparse.diags(x) @ matrix @ sparse.diags(x)

def ice_normalize(matrix, max_iter=50):
    n = matrix.shape[0]
    m = matrix.copy()

    for _ in range(max_iter):
        s = m.sum(axis=0)
        mask = (s == 0)
        s[mask] = 1
        v = s.mean()
        m = (m / s) * v
        m = (m / s[:, None]) * v

        if np.allclose(s, v):
            break

    return m

def normalize_matrix(matrix, method):
    if method == "raw":
        return matrix
    
    bin_count = matrix["bin1"].max() + 1
    rows = matrix["bin1"].values
    cols = matrix["bin2"].values
    data = matrix["count"].values
    sparse_matrix = sparse.csr_matrix((data, (rows, cols)), shape=(bin_count, bin_count))
    
    sparse_matrix = sparse_matrix + sparse_matrix.T - sparse.diags(sparse_matrix.diagonal())
    
    if method == "KR":
        normalized = kr_normalize(sparse_matrix)
    elif method == "ICE":
        normalized = ice_normalize(sparse_matrix)
    else:
        raise ValueError(f"Unknown normalization method: {method}")
    
    normalized = normalized.tocoo()
    return pd.DataFrame({
        "bin1": normalized.row,
        "bin2": normalized.col,
        "count": normalized.data
    })

def process_matrix(matrix, resolution):
    matrix = matrix[abs(matrix["bin1"] - matrix["bin2"]) > 2 * resolution]
    
    if matrix.empty or matrix["count"].max() == 0:
        return matrix

    counts = matrix["count"].values
    top_10_percent = heapq.nlargest(max(1, int(len(counts) * 0.1)), counts)
    threshold = min(top_10_percent)
    matrix.loc[matrix["count"] > threshold, "count"] = threshold
    return matrix

def extract_genes(gff_file, gene_list):
    if not gff_file or not gene_list:
        return defaultdict(list)
    genes = defaultdict(list)
    with open(gff_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            if fields[2] == "gene":
                gene_info = dict(item.split("=") for item in fields[8].split(";"))
                if gene_info.get("Name") in gene_list:
                    genes[fields[0]].append((int(fields[3]), int(fields[4]), gene_info["Name"]))
    return genes

def load_centromeres(centromere_file):
    if not centromere_file:
        return {}
    centromeres = {}
    with open(centromere_file, 'r') as f:
        for line in f:
            chrom, start, end = line.strip().split('\t')
            centromeres[chrom] = (int(start), int(end))
    return centromeres

def plot_heatmap(matrix, chrom, resolution, genes, centromere, max_value, output_pdf):
    bin_count = matrix["bin1"].max() // resolution + 1
    heatmap = np.zeros((bin_count, bin_count))
    for _, row in matrix.iterrows():
        i, j = row["bin1"] // resolution, row["bin2"] // resolution
        heatmap[i, j] = heatmap[j, i] = row["count"]

    fig, ax = plt.subplots(figsize=(10, 10))
    im = ax.imshow(heatmap, cmap="YlOrRd", vmin=0, vmax=max_value)
    plt.colorbar(im)

    ax.set_xticks(np.arange(0, bin_count, 5))
    ax.set_yticks(np.arange(0, bin_count, 5))
    ax.set_xticklabels(np.arange(0, bin_count, 5))
    ax.set_yticklabels(np.arange(0, bin_count, 5))

    for gene_start, gene_end, gene_name in genes.get(chrom, []):
        gene_bin = gene_start // resolution
        ax.plot([gene_bin, gene_bin], [bin_count, bin_count+0.5], color='red', linewidth=0.5)
        ax.plot([bin_count, bin_count+0.5], [gene_bin, gene_bin], color='red', linewidth=0.5)

    if chrom in centromere:
        centro_start, centro_end = centromere[chrom]
        centro_bin_start, centro_bin_end = centro_start // resolution, centro_end // resolution
        ax.plot([centro_bin_start, centro_bin_end], [bin_count, bin_count], color='gray', linewidth=2)
        ax.plot([bin_count, bin_count], [centro_bin_start, centro_bin_end], color='gray', linewidth=2)

    ax.set_title(f"Chromosome {chrom}")
    plt.tight_layout()
    output_pdf.savefig(fig)
    plt.close(fig)

def main():
    args = parse_arguments()
    resolution = args.resolution
    chrom_bins, chrom_sizes = load_chrom_bins(args.chrom_sizes, resolution)
    genes = extract_genes(args.gff, args.gene_list) if args.gff and args.gene_list else defaultdict(list)
    centromeres = load_centromeres(args.centromeres) if args.centromeres else {}

    max_values = defaultdict(float)

    for matrix_file in args.matrices:
        for chrom in chrom_sizes["chrom"]:
            matrix = load_matrix(matrix_file, chrom, chrom_bins, resolution)
            if not matrix.empty:
                matrix = normalize_matrix(matrix, args.normalization)
                matrix = process_matrix(matrix, resolution)
                max_values[chrom] = max(max_values[chrom], matrix["count"].max())

    with PdfPages(args.output) as pdf:
        for matrix_file in args.matrices:
            for chrom in chrom_sizes["chrom"]:
                matrix = load_matrix(matrix_file, chrom, chrom_bins, resolution)
                if not matrix.empty:
                    matrix = normalize_matrix(matrix, args.normalization)
                    matrix = process_matrix(matrix, resolution)
                    plot_heatmap(matrix, chrom, resolution, genes, centromeres, max_values[chrom], pdf)
                else:
                    print(f"Warning: No data for chromosome {chrom} in file {matrix_file}")

if __name__ == "__main__":
    main()
