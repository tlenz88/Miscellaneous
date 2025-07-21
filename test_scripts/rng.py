import random
from collections import Counter
import heapq
import os
import psutil
import tempfile
import argparse
import struct
import sys
import gc

def parse_arguments():
    parser = argparse.ArgumentParser(description="Efficient subsampling of large text files.")
    parser.add_argument("-i", "--input_file", help="Path to the input file")
    parser.add_argument("-n", "--num_reads", type=int, help="Number of reads to sample")
    parser.add_argument("-o", "--output", help="Path to the output file (default: input_file_subsample.ext)")
    parser.add_argument("-s", "--num_samples", type=int, default=100, help="Number of subsamples to generate (default: 100)")
    parser.add_argument("-t", "--temp_dir", help="Directory for temporary files (default: same as input file)")
    parser.add_argument("-m", "--memory", type=float, help="Fraction of available memory to use (default: 0.5)")
    parser.add_argument("-t", "--threads", type=int, default=0, help="Number of threads to use (0 for auto)")
    return parser.parse_args()

def get_default_output(input_file):
    base, ext = os.path.splitext(input_file)
    return f"{base}_subsample{ext}"

def estimate_memory_usage(num_reads, total_lines, num_samples):
    int_size = (total_lines.bit_length() + 7) // 8  # Round up to nearest byte
    return num_reads * int_size * num_samples

def get_available_memory():
    return psutil.virtual_memory().available


def run_rust_sampler(input_file, num_reads, num_samples, output_file, threads, mem):
    rs_script <- os.path.join(os.path.dirname(os.path.abspath(__file__)), "/target/release/avp_sampler")
    rust_command = [
        rs_script,
        input_file,
        "-n", str(num_reads),
        "-s", str(num_samples),
        "-o", output_file,
        "-t", str(threads)
    ]
    subprocess.run(rust_command, check=True)

def random_integer_sample_and_count(total_lines, num_reads, num_samples, max_memory, temp_dir):
    read_counts = Counter()
    estimated_memory = estimate_memory_usage(num_reads, total_lines, num_samples)
    
    # Determine the number of bytes needed to store each integer
    bytes_per_int = (total_lines.bit_length() + 7) // 8
    int_format = f'>I' if bytes_per_int <= 4 else '>Q'
    
    if estimated_memory > max_memory:
        print("Estimated memory usage exceeds limit. Using file-based approach.")
        temp_files = []
        try:
            for _ in range(num_samples):
                sample = random.sample(range(total_lines), num_reads)
                
                # Write each subsample to a temporary file
                with tempfile.NamedTemporaryFile(mode='wb', delete=False, dir=temp_dir) as temp_file:
                    temp_files.append(temp_file.name)
                    for num in sample:
                        temp_file.write(struct.pack(int_format, num))
                    temp_file.flush()  # Ensure data is flushed to disk
                    os.fsync(temp_file.fileno())  # Force the write to disk
                
                # Clear memory of the sample after writing to disk
                del sample
                gc.collect()  # Run garbage collection to free up memory
            
            # Process temporary files
            for temp_file_name in temp_files:
                with open(temp_file_name, 'rb') as f:
                    while True:
                        int_data = f.read(struct.calcsize(int_format))
                        if not int_data:
                            break
                        num = struct.unpack(int_format, int_data)[0]
                        read_counts[num] += 1
        finally:
            for temp_file_name in temp_files:
                os.unlink(temp_file_name)
    else:
        print("Sufficient memory available. Using in-memory approach.")
        for _ in range(num_samples):
            sample = random.sample(range(total_lines), num_reads)
            read_counts.update(sample)
    
    return read_counts

def create_consensus_sample(read_counts, num_reads, input_filename, output_filename):
    top_reads = set(heapq.nlargest(num_reads, read_counts, key=read_counts.get))
    
    with open(input_filename, 'r') as infile, open(output_filename, 'w') as outfile:
        for i, line in enumerate(infile):
            if i in top_reads:
                outfile.write(line)
                top_reads.remove(i)
                if not top_reads:
                    break

def count_lines(filename):
    with open(filename, 'rb') as f:
        return sum(1 for _ in f)

def main():
    args = parse_arguments()
    
    input_file = args.input_file
    num_reads = args.num_reads
    output_file = args.output or get_default_output(input_file)
    num_samples = args.num_samples
    temp_dir = args.temp_dir or os.path.dirname(input_file)
    memory_fraction = args.memory or 0.5
    
    available_memory = get_available_memory()
    max_memory = int(available_memory * memory_fraction)
    
    #total_lines = count_lines(input_file)
    total_lines = 388403417
    print(f"Total lines in input file: {total_lines}")
    
    #run_rust_sampler(args.input_file, args.num_reads, args.num_samples, rust_output, args.threads)

    #read_counts = random_integer_sample_and_count(total_lines, num_reads, num_samples, max_memory, temp_dir)
    create_consensus_sample(read_counts, num_reads, input_file, output_file)
    
    print(f"Consensus subsample written to: {output_file}")

if __name__ == '__main__':
    main()
