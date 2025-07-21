import sys
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages

# List of files and corresponding sample names
files = [sys.argv[1], sys.argv[2], sys.argv[3]]  # Add as many files as needed
sample_names = ['6hpi_vs_HFF_upreg', '24hpi_vs_6hpi_upreg', '24hpi_vs_6hpi_downreg']  # Corresponding sample names

# Load the list of pathways from File 2 (common pathways file)
pathway_file = 'immune_system_pathways_filtered.txt'
pathways_to_include = pd.read_csv(pathway_file, sep='\t', header=None)
pathway_list = pathways_to_include[0].tolist()  # First column is the list of pathways

# Initialize an empty DataFrame to store all results
combined_df = pd.DataFrame()

# Loop through each file and append the data
for i, file in enumerate(files):
    # Load ReactomePA output into a DataFrame for the current sample
    reactome_df = pd.read_csv(file, sep='\t')
    reactome_df['Sample'] = sample_names[i]  # Add a column for the sample name
    
    # Perform a left join to include all pathways from File 2, even if they aren't in the current sample
    merged_df = pd.DataFrame({'Description': pathway_list}).merge(
        reactome_df, on='Description', how='left'
    )
    
    # Convert qvalue to -log10(qvalue), handling missing values
    merged_df['-log10(qvalue)'] = -np.log10(merged_df['qvalue'].replace(0, np.nan))
    
    # Append the merged DataFrame to the combined DataFrame
    combined_df = pd.concat([combined_df, merged_df])

# Sort the DataFrame by pathway names (Description column)
combined_df = combined_df.sort_values(by='Description')

pdf = PdfPages('test.pdf')

# Create the dotplot
plt.figure(figsize=(12, 10))
dotplot = sns.scatterplot(
    data=combined_df,
    x='Count',
    y='Description',
    size='-log10(qvalue)',  # Dot size is determined by the -log10(qvalue)
    hue='Sample',  # Dot color is determined by the sample
    sizes=(50, 300),  # Scale dot size
    palette='Set2',  # Color palette for different samples
    legend='full'
)

# Customize the plot
plt.title('Dotplot of Pathways Across Samples')
plt.xlabel('Count')
plt.ylabel('Pathways')
plt.legend(bbox_to_anchor=(1, 1), loc='upper left', title='Sample')
plt.tight_layout()

# Save plot
pdf.savefig()
pdf.close()
