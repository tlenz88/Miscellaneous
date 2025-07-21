import h5py
import numpy as np
import matplotlib.pyplot as plt
import argparse
from pathlib import Path
from scipy.sparse import csr_matrix
import hdf5plugin

def validate_data(values, label):
    """
    Print diagnostic information about the data at various stages.
    """
    if len(values) == 0:
        print(f"{label} is empty!")
        return
        
    print(f"\n{label}:")
    print(f"  Mean: {np.mean(values):.4f}")
    print(f"  Median: {np.median(values):.4f}")
    print(f"  Min: {np.min(values):.4f}")
    print(f"  Max: {np.max(values):.4f}")
    print(f"  # positive values: {np.sum(values > 0)}")
    print(f"  # negative values: {np.sum(values < 0)}")
    print(f"  # zero values: {np.sum(values == 0)}")
    print(f"  Total values: {len(values)}")

def get_bin_size(h5_file):
    """
    Determine bin size from H5 file by looking at the difference between consecutive start positions.
    Returns bin size in base pairs.
    """
    with h5py.File(h5_file, 'r') as f:
        starts = f['intervals/start_list'][:]
        if len(starts) < 2:
            return None
        diff = np.diff(starts)
        bin_size = np.median(diff)
        return int(bin_size)

def prepare_distance_bins(min_dist, max_dist):
    """
    Create logarithmically spaced bins from minimum distance to maximum distance
    """
    n_decades = int(np.log10(max_dist) - np.log10(min_dist)) + 1
    n_bins = n_decades * 10
    return np.logspace(np.log10(min_dist), np.log10(max_dist), n_bins)

def bin_data(distances, values, bin_edges):
    """
    Bin the data using pre-defined bin edges and calculate mean values for each bin.
    Returns binned distances, values, and counts per bin for validation.
    """
    if len(distances) == 0 or len(values) == 0:
        return np.array([]), np.array([]), np.array([])
    
    # Find which bin each distance belongs to
    indices = np.digitize(distances, bin_edges)
    
    binned_values = []
    binned_distances = []
    counts_per_bin = []
    
    for i in range(1, len(bin_edges)):
        mask = indices == i
        if np.any(mask):
            bin_values = values[mask]
            binned_values.append(np.mean(bin_values))
            binned_distances.append(np.mean(distances[mask]))
            counts_per_bin.append(len(bin_values))
            
            # Print detailed bin statistics for debugging
            if len(bin_values) > 0:
                bin_start = bin_edges[i-1]
                bin_end = bin_edges[i]
                print(f"\nBin {i} ({bin_start:.0f}-{bin_end:.0f} bp):")
                print(f"  Values in bin: {len(bin_values)}")
                print(f"  Mean: {np.mean(bin_values):.4f}")
                print(f"  # positive: {np.sum(bin_values > 0)}")
                print(f"  # negative: {np.sum(bin_values < 0)}")
    
    return np.array(binned_distances), np.array(binned_values), np.array(counts_per_bin)

def load_and_plot_interactions(h5_files, bin_size=None, max_depth=None, skip_diagonal=False, output_file=None):
    """
    Load differential interactions from H5 files and create a line plot of distance vs log2FC.
    """
    plt.figure(figsize=(10, 6))
    colors = plt.cm.tab10(np.linspace(0, 1, len(h5_files)))
    
    if bin_size is None:
        bin_size = get_bin_size(h5_files[0])
        if bin_size is None:
            raise ValueError("Could not determine bin size from file. Please specify using --bin-size")
    
    if max_depth is None:
        max_depth = 100_000_000  # Default to 100 Mb
    
    bin_edges = prepare_distance_bins(bin_size, max_depth)
    all_log_values = []
    
    for file_idx, (file_path, color) in enumerate(zip(h5_files, colors)):
        print(f"\nProcessing file {file_idx + 1}: {file_path}")
        
        with h5py.File(file_path, 'r') as f:
            data = f['matrix/data'][:]
            indices = f['matrix/indices'][:]
            indptr = f['matrix/indptr'][:]
            shape = f['matrix/shape'][:]
            
            # Validate raw data
            validate_data(data, "Raw matrix data")
            
            matrix = csr_matrix((data, indices, indptr), shape=shape)
            starts = f['intervals/start_list'][:]
            
            rows, cols = matrix.nonzero()
            values = matrix.data
            distances = np.abs(starts[cols] - starts[rows])
            
            # Apply filters
            mask = np.ones(len(distances), dtype=bool)
            mask &= distances >= bin_size
            mask &= distances <= max_depth
            
            if skip_diagonal:
                mask &= rows != cols
            
            mask &= np.isfinite(values)
            
            if not np.any(mask):
                print(f"Warning: No valid data points found in {file_path}")
                continue
            
            filtered_distances = distances[mask]
            filtered_values = values[mask]
            
            # Validate filtered data
            validate_data(filtered_values, "Filtered values")
            
            # Bin the data and get counts per bin
            binned_distances, binned_values, counts = bin_data(filtered_distances, filtered_values, bin_edges)
            
            if len(binned_distances) == 0:
                print(f"Warning: No valid binned data for {file_path}")
                continue
            
            # Validate binned data
            validate_data(binned_values, "Binned values")
            
            # Plot both raw and binned data distribution
            plt.figure()
            plt.hist(filtered_values, bins=100, alpha=0.5)
            plt.title(f'Value distribution - {Path(file_path).stem}')
            plt.xlabel('log2FC')
            plt.ylabel('Count')
            if output_file:
                hist_file = f"{Path(output_file).stem}_hist_{file_idx}.png"
                plt.savefig(hist_file)
            plt.close()
            
            # Plot the binned data
            plt.figure()
            plt.plot(binned_distances,
                    binned_values,
                    color=color,
                    label=f"{Path(file_path).stem} (mean={np.mean(filtered_values):.4f})")
    
    if not all_log_values:
        raise ValueError("No valid data to plot")
    
    plt.xscale('log')
    plt.xlabel('Distance (bp)')
    plt.ylabel('log2 Fold Change')
    plt.axhline(y=0, color='black', linestyle='--', alpha=0.5)
    plt.grid(True, which='major', alpha=0.3)
    plt.grid(True, which='minor', alpha=0.1)
    
    def format_bp(x, pos):
        if x >= 1e6:
            return f'{x/1e6:.0f}Mb'
        elif x >= 1e3:
            return f'{x/1e3:.0f}kb'
        return f'{x:.0f}bp'
    
    plt.gca().xaxis.set_major_formatter(plt.FuncFormatter(format_bp))
    plt.xlim(bin_size, max_depth)
    
    max_abs_value = max(abs(min(all_log_values)), abs(max(all_log_values)))
    plt.ylim(-max_abs_value, max_abs_value)
    
    plt.legend()
    
    if output_file:
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
    else:
        plt.show()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Plot differential interactions from H5 files')
    parser.add_argument('h5_files', nargs='+', help='H5 files containing differential interactions')
    parser.add_argument('--bin-size', type=int, help='Bin size in base pairs (will try to determine from file if not provided)')
    parser.add_argument('--max-depth', type=int, help='Maximum interaction distance to plot (in base pairs)')
    parser.add_argument('--skip-diagonal', action='store_true', help='Skip diagonal interactions')
    parser.add_argument('--output', help='Output file path')
    
    args = parser.parse_args()
    
    load_and_plot_interactions(
        args.h5_files,
        bin_size=args.bin_size,
        max_depth=args.max_depth,
        skip_diagonal=args.skip_diagonal,
        output_file=args.output
    )
