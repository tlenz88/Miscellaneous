import sys
import argparse
import pysam
import pandas as pd
import re
import numpy as np
import matplotlib.pyplot as plt
from sklearn.mixture import GaussianMixture
from scipy.stats import norm
from scipy.optimize import fsolve
from pathlib import Path


def read_gff(gff_file):
    df = pd.read_csv(gff_file, sep="\t", header=None, comment="#")
    df = df[(df[2] == "protein_coding_gene") | (df[2] == "ncRNA_gene") | (df[2] == "gene")]
    df.reset_index(drop=True, inplace=True)
    
    def extract_gene_id(attr):
        match = re.search(r"(?:ID|Name|gene_id)=?([^;\"\n]+)", attr)
        return match.group(1) if match else "N/A"
    
    bed = pd.DataFrame({
        "chrom": df[0],
        "start": df[3] - 1,
        "end": df[4],
        "name": df[8].apply(extract_gene_id),
        "score": 0,
        "strand": df[6]
    })
    return bed

def get_window(row, upstream, downstream):
    if upstream is None or downstream is None:
        return pd.Series([row["start"], row["end"]])
    
    strand = row.get("strand", "+")
    if strand == "+":
        tss = row["start"]
        win_start = max(0, tss - upstream)
        win_end = tss + downstream
    elif strand == "-":
        tss = row["end"]
        win_start = max(0, tss - downstream)
        win_end = tss + upstream
    else:
        tss = row["start"]
        win_start = max(0, tss - upstream)
        win_end = tss + downstream
    return pd.Series([win_start, win_end])

def get_total_mapped_reads(bam_file):
    """Get total number of mapped reads from BAM file"""
    bam = pysam.AlignmentFile(bam_file, "rb")
    total_mapped = 0
    
    # Try to get from index stats first (faster)
    try:
        stats = bam.get_index_statistics()
        total_mapped = sum([stat.mapped for stat in stats])
    except:
        # Fallback: count manually (slower but more reliable)
        print(f"  Counting reads manually for {Path(bam_file).stem}...")
        for read in bam.fetch():
            if not read.is_unmapped:
                total_mapped += 1
    
    bam.close()
    return total_mapped

def calculate_gene_lengths(bed_df, upstream, downstream):
    """Calculate effective gene lengths for RPKM/FPKM normalization"""
    if upstream is not None and downstream is not None:
        # Using TSS windows - all have same length
        return np.full(len(bed_df), upstream + downstream)
    else:
        # Using gene body - calculate actual lengths
        return (bed_df["end"] - bed_df["start"]).values

def normalize_counts(count_df, sample_names, normalization_method, total_reads_dict, gene_lengths):
    """
    Normalize read counts using specified method
    
    Parameters:
    - count_df: DataFrame with raw counts
    - sample_names: List of sample names
    - normalization_method: Method to use for normalization
    - total_reads_dict: Dictionary mapping sample names to total read counts
    - gene_lengths: Array of gene/window lengths
    
    Returns:
    - DataFrame with normalized counts
    """
    normalized_df = count_df.copy()
    
    if normalization_method is None:
        return normalized_df
    
    print(f"\nApplying {normalization_method} normalization...")
    
    for sample_name in sample_names:
        count_col = f"{sample_name}_count"
        norm_col = f"{sample_name}_norm"
        
        raw_counts = count_df[count_col].values
        total_reads = total_reads_dict[sample_name]
        
        if normalization_method.upper() == "CPM":
            # Counts Per Million
            normalized_counts = (raw_counts / total_reads) * 1e6
            
        elif normalization_method.upper() == "RPKM":
            # Reads Per Kilobase per Million mapped reads
            gene_lengths_kb = gene_lengths / 1000.0
            rpk = raw_counts / gene_lengths_kb  # Reads per kilobase
            normalized_counts = rpk / (total_reads / 1e6)  # Per million mapped reads
            
        elif normalization_method.upper() == "FPKM":
            # Fragments Per Kilobase per Million (same as RPKM for single-end)
            gene_lengths_kb = gene_lengths / 1000.0
            fpk = raw_counts / gene_lengths_kb
            normalized_counts = fpk / (total_reads / 1e6)
            
        elif normalization_method.upper() == "TPM":
            # Transcripts Per Million
            gene_lengths_kb = gene_lengths / 1000.0
            rpk = raw_counts / gene_lengths_kb
            scaling_factor = np.sum(rpk) / 1e6
            normalized_counts = rpk / scaling_factor
            
        elif normalization_method.upper() == "QUANTILE":
            # Quantile normalization (rank-based)
            from scipy.stats import rankdata
            ranks = rankdata(raw_counts, method='average')
            normalized_counts = ranks / len(ranks) * np.max(raw_counts)
            
        elif normalization_method.upper() == "TMM":
            # Trimmed Mean of M-values (simplified version)
            # Calculate scaling factor relative to a reference sample
            if sample_name == sample_names[0]:
                # Use first sample as reference
                normalized_counts = raw_counts
            else:
                ref_counts = count_df[f"{sample_names[0]}_count"].values
                # Calculate M-values (log2 fold change)
                valid_mask = (raw_counts > 0) & (ref_counts > 0)
                if np.sum(valid_mask) > 0:
                    m_values = np.log2(raw_counts[valid_mask] / ref_counts[valid_mask])
                    # Trim extreme values (top and bottom 30%)
                    trimmed_mean = np.mean(np.sort(m_values)[int(0.3*len(m_values)):int(0.7*len(m_values))])
                    scaling_factor = 2 ** (-trimmed_mean)
                    normalized_counts = raw_counts * scaling_factor
                else:
                    normalized_counts = raw_counts
        
        else:
            raise ValueError(f"Unknown normalization method: {normalization_method}")
        
        # Add normalized counts to dataframe
        normalized_df[norm_col] = normalized_counts
        
        print(f"  {sample_name}: Total reads = {total_reads:,}, "
              f"Mean raw = {np.mean(raw_counts):.2f}, "
              f"Mean normalized = {np.mean(normalized_counts):.2f}")
    
    return normalized_df

def count_reads_single_bam(bam_file, bed_df, upstream, downstream):
    """Count reads for a single BAM file"""
    bam = pysam.AlignmentFile(bam_file, "rb")
    windows = bed_df.apply(get_window, axis=1, upstream=upstream, downstream=downstream)
    windows.columns = ["win_start", "win_end"]

    total = len(bed_df)
    read_counts = []
    sample_name = Path(bam_file).stem
    
    print(f"Processing {sample_name}...")
    for i, row in enumerate(bed_df.join(windows).itertuples(index=False), 1):
        count = bam.count(contig=row.chrom, start=int(row.win_start), end=int(row.win_end))
        read_counts.append(count)
        if i % 100 == 0 or i == total:
            print(f"  {sample_name}: {i} / {total} genes", end='\r', flush=True)
    print(f"  {sample_name}: {total} / {total} genes - Complete!")
    bam.close()
    return read_counts, sample_name

def count_reads_multiple_bams(bam_files, bed_df, upstream, downstream, normalization_method=None):
    """Count reads for multiple BAM files and return combined DataFrame with optional normalization"""
    # Start with the base gene information
    result_df = bed_df.copy()
    all_counts = []
    sample_names = []
    total_reads_dict = {}
    
    # Get total read counts for normalization if needed
    if normalization_method:
        print("\nCalculating total mapped reads for normalization...")
        for bam_file in bam_files:
            sample_name = Path(bam_file).stem
            total_reads = get_total_mapped_reads(bam_file)
            total_reads_dict[sample_name] = total_reads
            print(f"  {sample_name}: {total_reads:,} mapped reads")
    
    # Count reads per gene/region
    for bam_file in bam_files:
        read_counts, sample_name = count_reads_single_bam(bam_file, bed_df, upstream, downstream)
        result_df[f"{sample_name}_count"] = read_counts
        all_counts.extend(read_counts)
        sample_names.append(sample_name)
        if not normalization_method:  # Only add to dict if not already done
            total_reads_dict[sample_name] = len(read_counts)  # Placeholder
    
    # Apply normalization if requested
    if normalization_method:
        gene_lengths = calculate_gene_lengths(bed_df, upstream, downstream)
        result_df = normalize_counts(result_df, sample_names, normalization_method, 
                                   total_reads_dict, gene_lengths)
        
        # Update all_counts to use normalized values for threshold calculation
        all_counts = []
        for sample_name in sample_names:
            norm_col = f"{sample_name}_norm"
            all_counts.extend(result_df[norm_col].tolist())
    
    # Add combined statistics (using normalized counts if available)
    if normalization_method:
        count_columns = [f"{name}_norm" for name in sample_names]
        suffix = "_norm"
    else:
        count_columns = [f"{name}_count" for name in sample_names]
        suffix = "_count"
    
    result_df[f"total{suffix}"] = result_df[count_columns].sum(axis=1)
    result_df[f"mean{suffix}"] = result_df[count_columns].mean(axis=1)
    result_df[f"max{suffix}"] = result_df[count_columns].max(axis=1)
    
    return result_df, all_counts, sample_names

def calculate_threshold(all_counts, sample_names, plot_title="Multi-sample GMM Analysis", normalization_method=None, gmm_shift=0.0):
    """Calculate threshold using all counts from all samples"""
    counts = np.array(all_counts)
    counts = counts[counts >= 0]
    
    norm_text = f" ({normalization_method} normalized)" if normalization_method else " (raw counts)"
    print(f"\nAnalyzing {len(counts)} total count observations from {len(sample_names)} samples{norm_text}:")
    for name in sample_names:
        print(f"  - {name}")
    
    # Add pseudocount and log-transform
    log_counts = np.log1p(counts)

    # Reshape and fit GMM in log-space
    X_log = log_counts.reshape(-1, 1)
    gmm = GaussianMixture(n_components=2, random_state=0)
    gmm.fit(X_log)

    # Extract parameters
    means = gmm.means_.flatten()
    stds = np.sqrt(gmm.covariances_.flatten())
    weights = gmm.weights_.flatten()
    order = np.argsort(means)
    mean1, mean2 = means[order]
    std1, std2 = stds[order]
    weight1, weight2 = weights[order]

    # PDF functions in log-space
    def gauss1(x): return weight1 * norm.pdf(x, mean1, std1)
    def gauss2(x): return weight2 * norm.pdf(x, mean2, std2)
    def difference(x): return gauss1(x) - gauss2(x)

    # Find cutoff in log-space
    initial_guess = (mean1 + mean2) / 2
    cutoff_log = fsolve(difference, initial_guess)[0]
    
    # Apply shift in log-space (negative shift = lower threshold = more genes)
    adjusted_cutoff_log = cutoff_log + gmm_shift
    cutoff = np.expm1(cutoff_log)
    adjusted_cutoff = np.expm1(adjusted_cutoff_log)

    count_type = f"{normalization_method} normalized count" if normalization_method else "raw read count"
    
    print(f"\nGMM Analysis Results:")
    print(f"  Component 1 (Noise): μ={mean1:.3f}, σ={std1:.3f}, weight={weight1:.3f}")
    print(f"  Component 2 (Signal): μ={mean2:.3f}, σ={std2:.3f}, weight={weight2:.3f}")
    print(f"  Base GMM intersection ({count_type}): {cutoff:.2f}")
    if gmm_shift != 0.0:
        print(f"  Shift applied: {gmm_shift:.2f} (log-space)")
        print(f"  Adjusted cutoff ({count_type}): {adjusted_cutoff:.2f}")
    else:
        print(f"  Suggested cutoff ({count_type}): {cutoff:.2f}")

    # Create visualization
    fig, axes = plt.subplots(2, 2, figsize=(15, 12))
    
    # Plot 1: Combined log-scale histogram with GMM
    ax1 = axes[0, 0]
    x_vals = np.linspace(log_counts.min(), log_counts.max(), 1000)
    ax1.hist(log_counts, bins=100, density=True, alpha=0.4, label=f"Log({count_type})", color="gray")
    ax1.plot(x_vals, gauss1(x_vals), label="GMM Component 1 (Noise)", color="blue", linewidth=2)
    ax1.plot(x_vals, gauss2(x_vals), label="GMM Component 2 (Signal)", color="green", linewidth=2)
    ax1.plot(x_vals, gauss1(x_vals) + gauss2(x_vals), label="Combined GMM", color="red", linestyle="--")
    ax1.axvline(cutoff_log, color="orange", linestyle=":", label=f'Base intersection = {cutoff:.2f}', linewidth=2, alpha=0.7)
    if gmm_shift != 0.0:
        ax1.axvline(adjusted_cutoff_log, color="red", linestyle="--", label=f'Adjusted cutoff = {adjusted_cutoff:.2f}', linewidth=2)
    else:
        ax1.axvline(cutoff_log, color="red", linestyle="--", label=f'Cutoff = {cutoff:.2f}', linewidth=2)
    ax1.set_xlabel(f"Log({count_type} + 1)")
    ax1.set_ylabel("Density")
    ax1.set_title(f"{plot_title}\nLog-space GMM Fit")
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    
    # Plot 2: Linear scale histogram
    ax2 = axes[0, 1]
    ax2.hist(counts, bins=100, density=True, alpha=0.4, color="lightblue")
    threshold_to_show = adjusted_cutoff if gmm_shift != 0.0 else cutoff
    ax2.axvline(threshold_to_show, color="red", linestyle="--", label=f'Cutoff = {threshold_to_show:.2f}', linewidth=2)
    if gmm_shift != 0.0:
        ax2.axvline(cutoff, color="orange", linestyle=":", label=f'Base = {cutoff:.2f}', linewidth=2, alpha=0.7)
    ax2.set_xlabel(count_type.title())
    ax2.set_ylabel("Density")
    ax2.set_title("Linear Scale Distribution")
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    
    # Plot 3: Sample-wise statistics
    ax3 = axes[1, 0]
    n_genes = len(counts) // len(sample_names)
    sample_counts = [counts[i*n_genes:(i+1)*n_genes] for i in range(len(sample_names))]
    
    box_data = [sample_count for sample_count in sample_counts]
    bp = ax3.boxplot(box_data, tick_labels=sample_names, patch_artist=True)
    for patch in bp['boxes']:
        patch.set_facecolor('lightblue')
    threshold_to_show = adjusted_cutoff if gmm_shift != 0.0 else cutoff
    ax3.axhline(threshold_to_show, color="red", linestyle="--", label=f'Cutoff = {threshold_to_show:.2f}', linewidth=2)
    if gmm_shift != 0.0:
        ax3.axhline(cutoff, color="orange", linestyle=":", label=f'Base = {cutoff:.2f}', linewidth=2, alpha=0.7)
    ax3.set_ylabel(count_type.title())
    ax3.set_title(f"Per-sample {count_type.title()} Distribution")
    ax3.legend()
    ax3.grid(True, alpha=0.3)
    plt.setp(ax3.get_xticklabels(), rotation=45, ha='right')
    
    # Plot 4: Filtering impact
    ax4 = axes[1, 1]
    genes_kept_per_sample = []
    total_genes = n_genes
    
    for i, sample_count in enumerate(sample_counts):
        kept = np.sum(sample_count >= threshold_to_show)
        genes_kept_per_sample.append(kept)
    
    bars = ax4.bar(sample_names, genes_kept_per_sample, color='lightgreen', alpha=0.7)
    ax4.axhline(total_genes, color="gray", linestyle="--", label=f'Total genes = {total_genes}', alpha=0.7)
    ax4.set_ylabel("Genes Kept")
    ax4.set_title("Genes Retained After Filtering")
    ax4.legend()
    ax4.grid(True, alpha=0.3)
    plt.setp(ax4.get_xticklabels(), rotation=45, ha='right')
    
    # Add percentage labels on bars
    for bar, kept in zip(bars, genes_kept_per_sample):
        height = bar.get_height()
        percent = (kept / total_genes) * 100
        ax4.text(bar.get_x() + bar.get_width()/2., height + total_genes*0.01,
                f'{percent:.1f}%', ha='center', va='bottom', fontsize=9)
    
    plt.tight_layout()
    plt.show()
    
    return adjusted_cutoff if gmm_shift != 0.0 else cutoff

def filter_by_threshold(bed_df, threshold, filter_method="any", normalization_method=None):
    """
    Filter genes by threshold with different strategies:
    - "any": Keep if ANY sample meets threshold
    - "all": Keep if ALL samples meet threshold  
    - "mean": Keep if MEAN across samples meets threshold
    - "total": Keep if TOTAL across samples meets threshold
    """
    # Determine which columns to use based on normalization
    if normalization_method:
        count_suffix = "_norm"
        exclude_cols = {'total_norm', 'mean_norm', 'max_norm', 'total_count', 'mean_count', 'max_count'}
    else:
        count_suffix = "_count"
        exclude_cols = {'total_count', 'mean_count', 'max_count'}
    
    if filter_method == "any":
        # Keep if any sample has >= threshold reads
        count_cols = [col for col in bed_df.columns if col.endswith(count_suffix) and col not in exclude_cols]
        mask = (bed_df[count_cols] >= threshold).any(axis=1)
    elif filter_method == "all":
        # Keep if all samples have >= threshold reads
        count_cols = [col for col in bed_df.columns if col.endswith(count_suffix) and col not in exclude_cols]
        mask = (bed_df[count_cols] >= threshold).all(axis=1)
    elif filter_method == "mean":
        # Keep if mean across samples >= threshold
        mean_col = f"mean{count_suffix}"
        mask = bed_df[mean_col] >= threshold
    elif filter_method == "total":
        # Keep if total across samples >= threshold
        total_col = f"total{count_suffix}"
        mask = bed_df[total_col] >= threshold
    else:
        raise ValueError("filter_method must be one of: 'any', 'all', 'mean', 'total'")
    
    filtered_df = bed_df[mask].copy()
    norm_text = f" ({normalization_method})" if normalization_method else ""
    print(f"\nFiltering results ({filter_method} method{norm_text}):")
    print(f"  Threshold: {threshold:.2f}")
    print(f"  Genes kept: {len(filtered_df)} of {len(bed_df)} ({len(filtered_df)/len(bed_df)*100:.1f}%)")
    
    return filtered_df

def write_bed(filtered_df, output_file):
    """Write filtered genes to BED format"""
    filtered_df[["chrom", "start", "end", "name"]].to_csv(output_file, sep="\t", header=False, index=False)

def write_counts_table(filtered_df, output_file):
    """Write detailed count table"""
    filtered_df.to_csv(output_file, sep="\t", index=False)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Filter genes by read count across multiple BAM files using GMM with optional normalization.")
    parser.add_argument("--bam", nargs='+', required=True, help="Input BAM files (space-separated)")
    parser.add_argument("--gff", required=True, help="Input GFF file of genes")
    parser.add_argument("--threshold", type=float, help="Manual threshold (optional - will use GMM if not provided)")
    parser.add_argument("--upstream", type=int, help="Distance upstream of TSS (optional)")
    parser.add_argument("--downstream", type=int, help="Distance downstream of TSS (optional)")
    parser.add_argument("--normalize", choices=["CPM", "RPKM", "FPKM", "TPM", "QUANTILE", "TMM"], 
                        help="Normalization method: CPM (Counts Per Million), RPKM (Reads Per Kilobase Million), "
                             "FPKM (Fragments Per Kilobase Million), TPM (Transcripts Per Million), "
                             "QUANTILE (Quantile normalization), TMM (Trimmed Mean of M-values). "
                             "If not specified, raw counts are used.")
    parser.add_argument("--gmm-shift", type=float, default=0.0,
                        help="Shift GMM cutoff in log-space. Negative values lower the threshold (include more genes), "
                             "positive values raise it (exclude more genes). Example: -0.5 for more relaxed filtering.")
    parser.add_argument("--filter-method", choices=["any", "all", "mean", "total"], default="any",
                        help="Filtering strategy: 'any'=any sample ≥ threshold, 'all'=all samples ≥ threshold, 'mean'=mean ≥ threshold, 'total'=sum ≥ threshold")
    parser.add_argument("--out", default="filtered_genes.bed", help="Output BED file")
    parser.add_argument("--counts-out", help="Output detailed counts table (optional)")
    parser.add_argument("--auto-threshold", action="store_true", help="Calculate GMM threshold automatically and exit")
    
    args = parser.parse_args()

    # Validate BAM files
    for bam_file in args.bam:
        if not Path(bam_file).exists():
            print(f"Error: BAM file not found: {bam_file}")
            sys.exit(1)

    print(f"Processing {len(args.bam)} BAM files:")
    for bam_file in args.bam:
        print(f"  - {bam_file}")

    if args.normalize:
        print(f"Normalization method: {args.normalize}")

    # Read gene annotations
    bed = read_gff(args.gff)
    print(f"Loaded {len(bed)} genes from {args.gff}")

    # Count reads across all samples with optional normalization
    counted_df, all_counts, sample_names = count_reads_multiple_bams(
        args.bam, bed, args.upstream, args.downstream, args.normalize)

    # Calculate GMM threshold
    region_desc = (
        f"TSS -{args.upstream or 0}/+{args.downstream or 0}"
        if args.upstream and args.downstream
        else "gene body"
    )
    
    cutoff = calculate_threshold(all_counts, sample_names, 
                               plot_title=f"Multi-sample Analysis ({region_desc})",
                               normalization_method=args.normalize,
                               gmm_shift=args.gmm_shift or 0.0)

    # Use manual threshold if provided, otherwise use GMM
    final_threshold = args.threshold if args.threshold is not None else cutoff
    
    if args.auto_threshold:
        norm_text = f" ({args.normalize} normalized)" if args.normalize else ""
        print(f"\nRecommended threshold{norm_text}: {cutoff:.2f}")
        print("Use this value with --threshold parameter for actual filtering.")
        sys.exit(0)
    
    # Filter genes
    filtered = filter_by_threshold(counted_df, final_threshold, args.filter_method, args.normalize)
    
    # Write outputs
    write_bed(filtered, args.out)
    if args.counts_out:
        write_counts_table(filtered, args.counts_out)
    
    print(f"\nOutput files:")
    print(f"  BED file: {args.out}")
    if args.counts_out:
        print(f"  Counts table: {args.counts_out}")
    
    norm_text = f" ({args.normalize})" if args.normalize else ""
    print(f"\nSummary:")
    print(f"  Region analyzed: {region_desc}")
    print(f"  Samples: {len(sample_names)}")
    print(f"  Normalization: {args.normalize or 'None'}")
    print(f"  Filter method: {args.filter_method}")
    print(f"  Threshold used{norm_text}: {final_threshold:.2f}")
    print(f"  Genes retained: {len(filtered)} / {len(bed)} ({len(filtered)/len(bed)*100:.1f}%)")
