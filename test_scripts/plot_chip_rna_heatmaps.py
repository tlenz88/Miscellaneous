import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import os
import sys

# === INPUTS ===
chipseq_files = [
    "/mnt/f/ready_to_be_compressed/loic_project/ChIPseq/output_files/HFF_D3_with_genes.bed", 
    "/mnt/f/ready_to_be_compressed/loic_project/ChIPseq/output_files/HFF_D5_with_genes.bed", 
    "/mnt/f/ready_to_be_compressed/loic_project/ChIPseq/output_files/HFF_D7_with_genes.bed", 
    "/mnt/f/ready_to_be_compressed/loic_project/ChIPseq/output_files/Astro_D3_with_genes.bed", 
    "/mnt/f/ready_to_be_compressed/loic_project/ChIPseq/output_files/Astro_D5_with_genes.bed", 
    "/mnt/f/ready_to_be_compressed/loic_project/ChIPseq/output_files/Astro_D7_with_genes.bed"
]

# Initialize variables
gff_df = None
rna_df = None
gene_lengths = None

# === PROCESS COMMAND LINE ARGUMENTS ===
for file in sys.argv[1:]:
    ext = os.path.splitext(file)[1].lower()

    if ext == ".gff":
        print(f"Loading GFF file: {file}")
        gff_df = pd.read_csv(file, sep="\t", header=None, comment="#")
        gff_df.columns = ["Chromosome", "Source", "Feature", "Start", "End", "Score", "Strand", "Frame", "Attributes"]
        genes_df = gff_df[gff_df["Feature"].str.lower().isin(["gene", "protein_coding_gene", "ncrna_gene"])].copy()

        # Extract gene name
        def extract_gene_id(attr):
            fields = dict(item.split("=", 1) for item in attr.strip().split(";") if "=" in item)
            return fields.get("ID") or fields.get("gene_id") or fields.get("Name") or "NA"

        genes_df["Gene"] = genes_df["Attributes"].apply(extract_gene_id)
        genes_df = genes_df.drop(columns=["Attributes"])
        genes_df.dropna(subset=["Gene"], inplace=True)

        gene_lengths = genes_df["End"] - genes_df["Start"]
        gene_lengths.index = genes_df["Gene"]

    elif ext == ".txt":
        print(f"Loading RNA-seq file: {file}")
        # === LOAD RNA-SEQ DATA ===
        rna_df = pd.read_csv(file, sep="\t")
        rna_df.columns = rna_df.columns.str.strip()
        rna_df = rna_df.rename(columns={'Gene_ID': 'Gene'})
        rna_df = rna_df.set_index('Gene')

# === Check that required files were loaded ===
if gene_lengths is None:
    print("Error: No GFF file provided or processed")
    sys.exit(1)
    
if rna_df is None:
    print("Error: No RNA-seq file provided or processed")
    sys.exit(1)

# === Ensure matching gene lists ===
common_genes = rna_df.index.intersection(gene_lengths.index)
print(f"Found {len(common_genes)} common genes between RNA-seq and GFF")

rna_df = rna_df.loc[common_genes]
gene_lengths = gene_lengths.loc[common_genes]

# === Normalize to RPKM ===
total_reads = rna_df.sum(axis=0)
rpkm_df = (rna_df.T * 1e9 / gene_lengths).T / total_reads

# === Average replicates per condition ===
conditions = ["HFF_D3", "HFF_D5", "HFF_D7", "Astro_D3", "Astro_D5", "Astro_D7"]
rna_avg = pd.DataFrame(index=rpkm_df.index)

for cond in conditions:
    rep_cols = [col for col in rpkm_df.columns if col.startswith(cond)]
    if not rep_cols:
        print(f"Warning: no replicates found for {cond}")
        # Create empty column with zeros
        rna_avg[cond] = 0
    else:
        rna_avg[cond] = rpkm_df[rep_cols].mean(axis=1)

# === Output directory ===
output_dir = "heatmaps"
os.makedirs(output_dir, exist_ok=True)

# === PLOT RNA HEATMAPS MATCHED TO CHIP GENE ORDER ===
for chip_file in chipseq_files:
    if not os.path.exists(chip_file):
        print(f"Warning: ChIP-seq file not found: {chip_file}")
        continue
        
    condition_name = os.path.basename(chip_file).replace("_with_genes.bed", "")
    print(f"Processing {condition_name}...")

    try:
        # Load ChIP-seq BED file
        chip_df = pd.read_csv(chip_file, sep='\t', header=0)
        
        # Check if 'Gene' column exists
        if 'Gene' not in chip_df.columns:
            print(f"Error: 'Gene' column not found in {chip_file}")
            print(f"Available columns: {list(chip_df.columns)}")
            continue
            
        gene_order = chip_df["Gene"].tolist()
        print(f"Found {len(gene_order)} genes in ChIP-seq file")

        # Match RNA-seq condition name
        if condition_name not in rna_avg.columns:
            print(f"Warning: RNA-seq condition {condition_name} not found")
            print(f"Available conditions: {list(rna_avg.columns)}")
            continue

        # Extract RNA-seq expression for this condition, ordered by ChIP-seq
        matched_rna = rna_avg[condition_name].reindex(gene_order).fillna(0)
        
        print(f"Matched {(matched_rna > 0).sum()} genes with non-zero expression")

        # Normalize for visualization (log-scale for clarity)
        rna_matrix = np.log1p(matched_rna.values.reshape(-1, 1))

        # Plot heatmap
        plt.figure(figsize=(3, 12))
        sns.heatmap(rna_matrix, cmap='Reds', yticklabels=False, 
                   cbar_kws={'label': 'log1p(RPKM)'})
        plt.title(f"{condition_name} RNA-seq")
        plt.xlabel("")
        plt.ylabel(f"Genes (n={len(gene_order)})")
        plt.tight_layout()
        
        output_file = os.path.join(output_dir, f"{condition_name}_RNA_heatmap.png")
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.close()
        
        print(f"Saved: {output_file}")
        
    except Exception as e:
        print(f"Error processing {chip_file}: {str(e)}")

print("Heatmap generation complete!")
