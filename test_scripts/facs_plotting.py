import argparse
import os
import sys
import re
import numpy as np
import matplotlib.pyplot as plt
from flowkit import Sample
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.path import Path

def parse_args(args):
    parser = argparse.ArgumentParser(description='Check the help flag')
    parser.add_argument('-f', '--flow_data', dest='flow_data', help='REQUIRED: .fcs file containing flow data.', required=True)
    parser.add_argument('-a', '--channelA', dest='channelA', help='REQUIRED: First channel from flow data.', nargs='+', required=True)
    parser.add_argument('-b', '--channelB', dest='channelB', help='REQUIRED: Second channel from flow data.', nargs='+', required=True)
    parser.add_argument('-g', '--gates', dest='gates', help='REQUIRED: File containing gate vertices.', required=True)
    parser.add_argument('-o', '--output', dest='output', help='Output PDF path and file name.', required=False)
    return parser.parse_args()

def linear_scale(data, min_value, max_value, bin_resolution=256):
    scaled_data = (data - min_value) / (max_value - min_value) * bin_resolution
    return np.clip(scaled_data, 0, bin_resolution)

def log_scale(data, min_value, max_value, bin_resolution=256):
    data = np.clip(data, 1e-10, None)
    scaled_data = np.log1p(data - min_value) / np.log1p(max_value - min_value) * bin_resolution
    return np.clip(scaled_data, 0, bin_resolution)

def sampling_and_scaling(sample, sample_channel):
    idx = sample.get_channel_index(sample_channel)
    data = sample.get_channel_events(idx, source='raw')
    pattern = re.compile(r'log', re.IGNORECASE)
    if pattern.search(sample_channel):
        return log_scale(data, data.min(), data.max())
    else:
        return linear_scale(data, data.min(), data.max())

def apply_polygonal_gate(sampleA, sampleB, vertices):
    points = np.column_stack((sampleA, sampleB))
    polygon_path = Path(vertices)
    mask = polygon_path.contains_points(points)
    return mask

def plot_channel(channelA, channelB, sampleA, sampleB, ax, vertices=None):
    hb = ax.hexbin(sampleA, sampleB, bins='log', gridsize=256, cmap='turbo_r', mincnt=1)
    cb = plt.colorbar(hb, ax=ax)
    cb.set_label('Density')

    ax.set_title(f"{channelA} vs {channelB}")
    ax.set_xlabel(channelA)
    ax.set_ylabel(channelB)
    
    ax.set_aspect('equal', adjustable='box')
    
    ax.set_xlim(0, 256)
    ax.set_ylim(0, 256)
    ticks = [0, 64, 128, 192, 256]
    ax.set_xticks(ticks)
    ax.set_yticks(ticks)

    if vertices is not None:
        ax.add_patch(plt.Polygon(vertices, fill=False, edgecolor='black'))

def read_gates(gates_file):
    gates = []
    current_gate = []
    with open(gates_file, 'r') as f:
        for line in f:
            line = line.strip()
            if line:
                if line.startswith('Gate'):
                    if current_gate:
                        gates.append(current_gate)
                        current_gate = []
                else:
                    x, y = map(float, line.split(','))
                    current_gate.append((x, y))
    if current_gate:
        gates.append(current_gate)
    return gates

def main():
    args = parse_args(sys.argv[1:])

    sample = Sample(args.flow_data)

    if args.output:
        pdf = PdfPages(args.output)
    else:
        output = os.path.splitext(os.path.abspath(args.flow_data))[0] + '.pdf'
        pdf = PdfPages(output)

    num_samples = len(args.channelA)
    fig, axes = plt.subplots(nrows=1, ncols=num_samples, figsize=(num_samples * 5, 5))

    if num_samples == 1:
        axes = [axes]

    gates = read_gates(args.gates)
    
    # Ensure we have enough gates for all channel pairs
    if len(gates) < num_samples:
        print(f"Warning: Not enough gates provided. Expected {num_samples}, got {len(gates)}.")
        gates.extend([None] * (num_samples - len(gates)))

    # Apply gates and plot sequentially
    mask = np.ones(len(sampling_and_scaling(sample, args.channelA[0])), dtype=bool)
    for i in range(num_samples):
        sampleA = sampling_and_scaling(sample, args.channelA[i])
        sampleB = sampling_and_scaling(sample, args.channelB[i])

        # Plot data (gated from previous iterations if i > 0)
        plot_channel(args.channelA[i], args.channelB[i], sampleA[mask], sampleB[mask], axes[i], gates[i])

        # Apply gate for next iteration if available
        if i < num_samples - 1 and gates[i] is not None:
            new_mask = apply_polygonal_gate(sampleA, sampleB, gates[i])
            mask = mask & new_mask

    plt.tight_layout()
    pdf.savefig()
    pdf.close()

if __name__ == '__main__':
    main()
