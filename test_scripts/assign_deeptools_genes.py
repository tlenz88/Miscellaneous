import pandas as pd
import pyranges as pr
import sys
import os
import re

for file in sys.argv[1:]:
    ext = os.path.splitext(file)[1].lower()

    if ext == ".bed":
        print(f"Loading BED file: {file}")
        bed_filename = file
        bed_df = pd.read_csv(file, sep="\t", header=None, comment="#")
        bed_df.columns = ["Chromosome", "Start", "End", "Name", "Score", "Strand", "ThickStart", "ThickEnd", "ItemRGB", "BlockCount", "BlockSizes", "BlockStart", "DeepTools_group"]
        bed_pr = pr.PyRanges(bed_df)

    elif ext == ".gtf":
        print(f"Loading GTF file: {file}")
        genes_df = pd.read_csv(file, sep="\t", header=None, comment="#")
        genes_df.columns = ["Chromosome", "Source", "Feature", "Start", "End", "Score", "Strand", "Frame", "Attributes"]
        def extract_gene_id(attr):
            match = re.search(r'gene_id\s+"([^"]+)"', attr)
            return match.group(1) if match else "NA"
        genes_df["Gene"] = genes_df["Attributes"].apply(extract_gene_id)
        genes_df = genes_df.drop(columns=["Attributes"])
        genes_pr = pr.PyRanges(genes_df)

overlap = bed_pr.join(genes_pr, apply_strand_suffix=False)
result_df = overlap.df
result_df["overlap_len"] = result_df.apply(lambda row: min(row["End"], row["End_b"]) - max(row["Start"], row["Start_b"]), axis=1)
cols_to_group = ["Chromosome", "Start", "End", "Name", "Score", "Strand"]
best_overlap_df = result_df.sort_values("overlap_len", ascending=False).drop_duplicates(subset=cols_to_group)
cols_to_keep = ["Chromosome", "Start", "End", "Name", "Score", "Strand_b", "ThickStart", "ThickEnd", "ItemRGB", "BlockCount", "BlockSizes", "BlockStart", "DeepTools_group", "Gene"]
final_df = best_overlap_df[cols_to_keep]
final_df.columns = ["Chromosome", "Start", "End", "Name", "Score", "Strand", "ThickStart", "ThickEnd", "ItemRGB", "BlockCount", "BlockSizes", "BlockStart", "DeepTools_group", "Gene"]

out_file = os.path.splitext(bed_filename)[0] + "_with_genes.bed"
final_df.to_csv(out_file, sep="\t", index=False)
