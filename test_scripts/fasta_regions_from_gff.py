#!/usr/bin/env python3

import argparse
from Bio import SeqIO
import pandas as pd

def parse_gff(gff_file):
    """Parse GFF and return a list of regions."""
    regions = []
    with open(gff_file) as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 9:
                continue
            chrom, source, feature_type, start, end, score, strand, phase, attributes = parts
            region_id = attributes.split(";")[0].split("=")[-1]
            regions.append({
                "chrom": chrom,
                "start": int(start),
                "end": int(end),
                "strand": strand,
                "id": region_id
            })
    return regions

def extract_sequences(fasta_file, regions):
    """Extract sequences from FASTA based on GFF regions."""
    fasta = SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta"))
    sequences = []
    for region in regions:
        chrom = region["chrom"]
        if chrom not in fasta:
            print(f"Warning: {chrom} not found in FASTA.")
            continue
        seq_record = fasta[chrom]
        subseq = seq_record.seq[region["start"] - 1 : region["end"]]
        if region["strand"] == "-":
            subseq = subseq.reverse_complement()
        header = f">{region['id']} {chrom}:{region['start']}-{region['end']}({region['strand']})"
        sequences.append((header, str(subseq)))
    return sequences

def write_fasta(sequences, output_file):
    """Write extracted sequences to a FASTA file."""
    with open(output_file, "w") as out:
        for header, seq in sequences:
            out.write(f"{header}\n")
            for i in range(0, len(seq), 60):
                out.write(seq[i:i+60] + "\n")

def main():
    parser = argparse.ArgumentParser(description="Extract sequences from FASTA using GFF regions.")
    parser.add_argument("-f", "--fasta", required=True, help="Input FASTA file")
    parser.add_argument("-g", "--gff", required=True, help="Input GFF file")
    parser.add_argument("-o", "--output", required=True, help="Output FASTA file")
    args = parser.parse_args()

    regions = parse_gff(args.gff)
    sequences = extract_sequences(args.fasta, regions)
    write_fasta(sequences, args.output)

if __name__ == "__main__":
    main()
