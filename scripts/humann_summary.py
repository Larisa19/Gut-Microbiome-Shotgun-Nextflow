
from pathlib import Path
import sys

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


# Input and output
input_file = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(
    "pathabundance_merged.tsv"
)
output_dir = Path(".")
output_dir.mkdir(exist_ok=True)

df = pd.read_csv(input_file, sep="\t")

pathway_column = "# Pathway"
if pathway_column not in df.columns:
    raise ValueError(f"Missing required column: {pathway_column}")

sample_columns = [
    col for col in df.columns if col.endswith("_Abundance")
]
if not sample_columns:
    raise ValueError("No sample abundance columns were found.")

# Retain identified pathways and exclude unmapped/unintegrated entries
pathways = df[
    ~df[pathway_column].str.startswith(
        ("UNMAPPED", "UNINTEGRATED"), na=False
    )
].copy()

pathways[sample_columns] = pathways[sample_columns].apply(
    pd.to_numeric, errors="raise"
)

# Save the filtered pathway table
pathways.to_csv(
    output_dir / "pathabundance_summary.tsv",
    sep="\t",
    index=False
)

# Select the 15 pathways with the highest mean abundance
pathways["mean_abundance"] = pathways[sample_columns].mean(axis=1)

top15 = pathways.nlargest(15, "mean_abundance").sort_values(
    "mean_abundance"
)

# Figure 1: overall ranking
fig, ax = plt.subplots(figsize=(10, 7))
ax.barh(top15[pathway_column], top15["mean_abundance"])
ax.set_xlabel("Mean pathway abundance across samples")
ax.set_ylabel("HUMAnN pathway")
ax.set_title("Top 15 metabolic pathways by mean abundance")
fig.tight_layout()

fig.savefig(
    output_dir / "top15_pathways.png",
    dpi=300,
    bbox_inches="tight"
)
plt.close(fig)

# Figure 2: heatmap across individual samples
# Keep the sample order consistent with the study design.
preferred_order = [
    "FBI016_W0_Abundance",
    "FBI016_W16_Abundance",
    "FBI029_W0_Abundance",
    "FBI029_W16_Abundance",
]
ordered_columns = [
    col for col in preferred_order if col in sample_columns
]
ordered_columns += [
    col for col in sample_columns if col not in ordered_columns
]

values = top15[ordered_columns].to_numpy(dtype=float)

# Log transform to make differences easier to visualize.
# Use half the smallest positive value as a pseudocount.
positive_values = values[values > 0]
pseudocount = (
    positive_values.min() / 2 if positive_values.size else 1.0
)
log_values = np.log10(values + pseudocount)

sample_labels = [
    col.removesuffix("_Abundance").replace("_", "\n")
    for col in ordered_columns
]

fig, ax = plt.subplots(figsize=(9, 8))
image = ax.imshow(
    log_values,
    aspect="auto",
    interpolation="nearest",
    cmap="viridis"
)

ax.set_xticks(range(len(sample_labels)))
ax.set_xticklabels(sample_labels)
ax.set_yticks(range(len(top15)))
ax.set_yticklabels(top15[pathway_column], fontsize=7)

ax.set_xlabel("Participant and timepoint")
ax.set_ylabel("HUMAnN pathway")
ax.set_title("Top 15 metabolic pathways across samples")

colorbar = fig.colorbar(image, ax=ax)
colorbar.set_label("log10(abundance + pseudocount)")

# Separate the two participants visually when all four samples exist.
if len(ordered_columns) == 4:
    ax.axvline(1.5, color="white", linewidth=2)

fig.tight_layout()
fig.savefig(
    output_dir / "top15_pathways_heatmap.png",
    dpi=300,
    bbox_inches="tight"
)
plt.close(fig)

print(f"Input pathways: {len(df)}")
print(f"Identified pathways retained: {len(pathways)}")
print(f"Top pathways selected: {len(top15)}")
print(f"Outputs written to: {output_dir}")
