#!/usr/bin/python3

import sys
import os
import random
import re
from Bio import SeqIO

def main():
    allValidPairs = []
    for arg in sys.argv[1:]:
        try:
            read_count = int(arg)
        except ValueError:
            allValidPairs.append(arg)
    try:
        print("Subsampling input .allValidPairs files to %s reads" % read_count)
    except NameError:
        print("An integer is required to filter reads.")
    sys.exit()
    for aVP in allValidPairs:
        with open(''.join([os.path.splitext(aVP)[0], '_subsample.allValidPairs']), "w") as infile:
            all_reads = list(SeqIO.parse(infile, 'fastq'))
            sampled_reads = random.sample(range(len(all_reads)), read_count)
        for idx in sampled_indices:
            SeqIO.write(sampled_reads[idx], R1, 'fastq')

if __name__ == "__main__":
    main()
