import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import ttest_ind
from statsmodels.stats.multitest import multipletests

# Load data
df = pd.read_csv("D7_replicates_profile.txt", sep="\t")
df.set_index("Sample", inplace=True)

# Infer condition from sample name
df["Condition"] = df.index.to_series().astype(str).apply(
    lambda x: "Astro" if "Astro" in x else ("HFF" if "HFF" in x else "Unknown")
)

# Convert Bin names to numeric
bin_columns = [col for col in df.columns if col != "Condition"]
rename_dict = {col: int(col.replace("Bin", "")) for col in bin_columns}
df.rename(columns=rename_dict, inplace=True)
bin_columns = sorted(rename_dict.values())

# Separate replicates
astro_df = df[df["Condition"] == "Astro"]
hff_df = df[df["Condition"] == "HFF"]

# Per-bin t-tests
pvals, log2fcs = [], []
for bin_id in bin_columns:
    a_vals = astro_df[bin_id].values
    h_vals = hff_df[bin_id].values
    if len(a_vals) >= 2 and len(h_vals) >= 2:
        _, p = ttest_ind(a_vals, h_vals, equal_var=False)
    else:
        p = np.nan
    log2fc = np.log2((np.mean(a_vals) + 1e-6) / (np.mean(h_vals) + 1e-6))
    pvals.append(p)
    log2fcs.append(log2fc)

# Multiple testing correction
valid = ~np.isnan(pvals)
adj_pvals = np.full_like(pvals, np.nan, dtype=np.float64)
if valid.sum() > 0:
    _, adj_valid, _, _ = multipletests(np.array(pvals)[valid], method="fdr_bh")
    adj_pvals = np.where(valid, adj_valid, np.nan)

# Save table
results = pd.DataFrame({
    "bin": bin_columns,
    "log2FC": log2fcs,
    "pval": pvals,
    "adj_pval": adj_pvals,
    "significant": adj_pvals < 0.05
})
results.to_csv("significant_bins.tsv", sep="\t", index=False)

# Plot
x = np.array(bin_columns)
plt.figure(figsize=(12, 6))
plt.plot(x, astro_df[bin_columns].mean(), label="Astro", color="blue")
plt.plot(x, hff_df[bin_columns].mean(), label="HFF", color="red")
for i, sig in enumerate(results["significant"]):
    if sig:
        plt.axvspan(x[i] - 0.5, x[i] + 0.5, color="green", alpha=0.3)
plt.xlabel("Bin (relative to TSS)")
plt.ylabel("Mean Signal")
plt.title("Astro vs HFF TSS Profile with Significant Bins Highlighted")
plt.legend()
plt.tight_layout()
plt.savefig("profile_with_significant_bins.png")
plt.close()
