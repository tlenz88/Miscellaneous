import gzip
import sys
import matplotlib.pyplot as plt
import numpy as np

def read_fastq(gzipped_fastq):
    """Read a gzipped FASTQ file and return a list of fragment lengths."""
    with gzip.open(gzipped_fastq, 'rt') as f:
        fragment_lengths = []
        for line in f:
            # Read four lines at a time (one entry in FASTQ)
            seq_line = next(f)  # Sequence line
            next(f)             # Plus line
            next(f)             # Quality line
            fragment_lengths.append(len(seq_line.strip()))
    return fragment_lengths

def calculate_statistics(fragment_lengths):
    """Calculate the statistics for fragment lengths."""
    total_reads = len(fragment_lengths)
    greater_than_50 = sum(length > 50 for length in fragment_lengths)
    greater_than_75 = sum(length > 75 for length in fragment_lengths)
    greater_than_100 = sum(length > 100 for length in fragment_lengths)
    greater_than_125 = sum(length > 125 for length in fragment_lengths)

    # Calculate percentages
    pct_greater_than_50 = (greater_than_50 / total_reads) * 100
    pct_greater_than_75 = (greater_than_75 / total_reads) * 100
    pct_greater_than_100 = (greater_than_100 / total_reads) * 100
    pct_greater_than_125 = (greater_than_125 / total_reads) * 100

    return (pct_greater_than_50, pct_greater_than_75, pct_greater_than_100, pct_greater_than_125)

def plot_fragment_distribution(fragment_lengths, stats, output_pdf):
    """Plot the distribution of fragment lengths and display statistics."""
    plt.figure(figsize=(10, 6))
    plt.hist(fragment_lengths, bins=range(0, max(fragment_lengths) + 1, 1), alpha=0.75, color='blue')
    
    plt.title('Fragment Size Distribution')
    plt.xlabel('Fragment Length (bp)')
    plt.ylabel('Number of Reads')
    plt.xlim(0, max(fragment_lengths) + 10)
    plt.grid(axis='y', alpha=0.75)

    # Adding percentage annotations to the plot
    plt.text(0.7 * max(fragment_lengths), 0.9 * plt.ylim()[1], f'> 50 bp: {stats[0]:.2f}%', fontsize=12, color='black')
    plt.text(0.7 * max(fragment_lengths), 0.85 * plt.ylim()[1], f'> 75 bp: {stats[1]:.2f}%', fontsize=12, color='black')
    plt.text(0.7 * max(fragment_lengths), 0.80 * plt.ylim()[1], f'> 100 bp: {stats[2]:.2f}%', fontsize=12, color='black')
    plt.text(0.7 * max(fragment_lengths), 0.75 * plt.ylim()[1], f'> 125 bp: {stats[3]:.2f}%', fontsize=12, color='black')

    # Save the plot to a PDF file
    plt.savefig(output_pdf, format='pdf')
    plt.close()  # Close the figure

def main(gzipped_fastq):
    fragment_lengths = read_fastq(gzipped_fastq)
    (pct_50, pct_75, pct_100, pct_125) = calculate_statistics(fragment_lengths)

    output_pdf = 'fragment_size_distribution.pdf'  # Define output PDF file name
    plot_fragment_distribution(fragment_lengths, (pct_50, pct_75, pct_100, pct_125), output_pdf)

# Example usage
if __name__ == "__main__":
    gzipped_fastq_file = sys.argv[1]  # Update with your file path
    main(gzipped_fastq_file)
