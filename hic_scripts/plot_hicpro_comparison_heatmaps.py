#!/usr/bin/env python3

import sys
import csv
import argparse
import numpy as np
import pandas as pd
import heapq
import re
import math
from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.colors import Normalize, LogNorm
import os


def parse_args(args):
    parser = argparse.ArgumentParser(description='Check the help flag')
    parser.add_argument('-ch',
                        '--chr_sizes',
                        dest='sizes',
                        help='REQUIRED: chromosome sizes (.chrom.sizes)',
                        required=True)
    parser.add_argument('-m1',
                        '--matrix1',
                        dest='matrix1',
                        help='REQUIRED: tab-delimited HiC-Pro matrix.',
                        required=True)
    parser.add_argument('-m2',
                        '--matrix2',
                        dest='matrix2',
                        help='REQUIRED: tab-delimited HiC-Pro matrix.',
                        required=True)
    parser.add_argument('-o',
                        '--output',
                        dest='output',
                        help='Output directory.',
                        required=True)
    parser.add_argument('-r',
                        '--resolution',
                        dest='res',
                        help='Binning resolution.',
                        required=False,
                        default=100000)
    parser.add_argument('-w',
                        '--window',
                        dest='window',
                        help='Size of sliding window to plot.',
                        required=False)
    return parser.parse_args()


def input_params(args):
    out = args.output
    if args.res is not None:
        res = resolution(args.res)
    if args.window is not None:
        window_size = int(args.window)
    else:
        window_size = res * 100
    sample = f'{re.match(rf"(.*?){str(res)}", 
    str(os.path.basename(args.matrix2))).group(1)}vs_{re.match(rf"(.*?){str(res)}", 
    str(os.path.basename(args.matrix1))).group(1)}'.strip('_')
    return out, sample, res, window_size


def resolution(res):
    res = res.lower()
    if 'kb' in res:
        return int(res.split('kb')[0]) * 1000
    elif 'mb' in res:
        return int(res.split('mb')[0]) * 1000000
    else:
        return int(res)


def get_chrom_starts(sizes, res):
    sizes[1] = sizes[1] // res + 1
    start_index = 0
    chrom_starts = {}
    for chrom in sizes.itertuples():
        chrom_starts[chrom[1]] = start_index
        start_index += chrom[2]
    return chrom_starts


def plot_intra(pdf,mat,maxColor,sample,ch,start,end,res):
    fig, ax = plt.subplots()
    #plt.imshow(mat, cmap='Reds', origin='lower', vmax=maxColor)
    plt.imshow(np.log1p(mat),
                cmap='Reds',
                origin='lower',
                norm=LogNorm())
    cbar = plt.colorbar(cmap='Reds')
    cbar.ax.set_ylabel('Normalized contact counts',
                        rotation=270, labelpad=12, size=6)
    plt.xticks(ticks=np.linspace(0, mat.shape[0] - 1, 6),
               labels=np.linspace(start, end, 6, dtype=float))
    plt.yticks(ticks=np.linspace(0, mat.shape[0] - 1, 6),
               labels=np.linspace(start, end, 6, dtype=float))
    plt.title(f'sample {ch}:{start}-{end} {int(res // 1000)}kb \n intrachromosomal interactions', 
                pad=12)
    plt.tight_layout()
    pdf.savefig()
    plt.close()


def main():
    args = parse_args(sys.argv[1:])
    out, sample, res, window_size = input_params(args)
    sizes = pd.read_csv(args.sizes, sep='\t', header=None)
    num_bins = np.sum(sizes[1] // res + 1)
    chrom_starts = get_chrom_starts(sizes, res)
    window_bins = window_size // res
    """
    mat = np.zeros((num_bins, num_bins))
    
    with open(args.matrix1, 'r') as m:
        for i in csv.reader(m, delimiter = '\t'):
            if int(i[0]) <= num_bins and int(i[1]) <= num_bins:
                mat[int(i[0]) - 1][int(i[1]) - 1] = i[2]

    with open(args.matrix2, 'r') as m:
        for i in csv.reader(m, delimiter = '\t'):
            if int(i[0]) <= num_bins and int(i[1]) <= num_bins:
                mat[int(i[1]) - 1][int(i[0]) - 1] = i[2]

    np.fill_diagonal(mat, 0)
    np.fill_diagonal(mat[:, 1:], 0)
    np.fill_diagonal(mat[1:, :], 0)
    """
    pdf = PdfPages(f'{out}/{sample}_{int(res / 1000)}kb_log1p.pdf')

    for c in chrom_starts:
        if c == 'chr1':
            continue
        chrom_length = sizes.loc[sizes[0] == c, 1].iloc[0]
        window_num = chrom_length // window_bins * 2 + 1

        print(f'Plotting chromosome: {c}')
        pdf = PdfPages(f'{out}/{sample}_{c}_{int(res / 1000)}kb_log1p.pdf')
        mat = np.zeros((chrom_length, chrom_length))

        with open(args.matrix1, 'r') as m:
            for line in csv.reader(m, delimiter = '\t'):
                if int(line[0]) > int(chrom_starts[c]) + chrom_length:
                    break
                elif int(line[0]) > int(chrom_starts[c]) and int(line[0]) <= int(chrom_starts[c]) + chrom_length and int(line[1]) > int(chrom_starts[c]) and int(line[1]) <=  int(chrom_starts[c]) + chrom_length:
                    mat[int(line[0]) - int(chrom_starts[c]) - 1][int(line[1]) - int(chrom_starts[c]) - 1] = line[2]

        with open(args.matrix2, 'r') as m:
            for line in csv.reader(m, delimiter = '\t'):
                if int(line[0]) > int(chrom_starts[c]) + chrom_length:
                    break
                elif int(line[0]) > int(chrom_starts[c]) and int(line[0]) <= int(chrom_starts[c]) + chrom_length and int(line[1]) > int(chrom_starts[c]) and int(line[1]) <=  int(chrom_starts[c]) + chrom_length:
                    mat[int(line[1]) - int(chrom_starts[c]) - 1][int(line[0]) - int(chrom_starts[c]) - 1] = line[2]

        np.fill_diagonal(mat, 0)
        np.fill_diagonal(mat[:, 1:], 0)
        np.fill_diagonal(mat[1:, :], 0)

        for w in range(window_num):
            #start = int(chrom_starts[c] + window_bins * w * 0.5)
            start = int(window_bins * w * 0.5)
            
            end = int(start + window_bins)
            #if end > chrom_starts[c] + chrom_length:
            #    end = chrom_starts[c] + chrom_length
            if end > chrom_length:
                end = chrom_length

            window_mat = mat[start:end, start:end]
            if np.sum(window_mat) == 0:
                continue

            maxColor = np.max(window_mat)
            plot_intra(pdf,
                       window_mat,
                       maxColor,
                       sample,
                       str(c),
                       start*res/1e6,#(start-chrom_starts[c])*res/1e6,
                       end*res/1e6,#(end-chrom_starts[c])*res/1e6,
                       res)
        pdf.close()
    #pdf.close()

if __name__ == '__main__':
    main()
