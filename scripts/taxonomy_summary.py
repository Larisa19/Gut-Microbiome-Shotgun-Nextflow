
import argparse
from pathlib import Path

import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def main():
    parser = argparse.ArgumentParser(
        description="Summarize MetaPhlAn genus-level abundances."
    )
    parser.add_argument("--input", required=True)
    parser.add_argument("--outdir", default=".")
    parser.add_argument("--top-n", type=int, default=15)
    args = parser.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    # Skip the MetaPhlAn database-version comment.
    df = pd.read_csv(args.input, sep="\t", skiprows=1)

    if "clade_name" not in df.columns:
        raise ValueError("Missing expected column: clade_name")

    sample_cols = [c for c in df.columns if c != "clade_name"]
    if not sample_cols:
        raise ValueError("No sample columns found.")

    # Select genus-level rows only, not species or SGB descendants.
    mask = df["clade_name"].astype(str).str.contains(
        r"(?:^|\|)g__[^|]+$", regex=True
    )
    genus = df.loc[mask].copy()

    genus["Genus"] = (
        genus["clade_name"]
        .str.split("|")
        .str[-1]
        .str.replace(r"^g__", "", regex=True)
    )
    genus = genus[~genus["Genus"].isin(["", "unclassified"])]

    if genus.empty:
        raise ValueError("No genus-level rows found.")

    abundance = genus.set_index("Genus")[sample_cols].apply(
        pd.to_numeric, errors="raise"
    )

    if (abundance < 0).any().any():
        raise ValueError("Negative relative abundances detected.")

    # Save all genus-level abundances as percentages.
    abundance.to_csv(
        outdir / "taxonomy_summary.tsv",
        sep="\t",
        float_format="%.6f"
    )

    # Rank genera by mean abundance across the selected samples.
    top = abundance.mean(axis=1).nlargest(args.top_n).index
    plot_data = abundance.loc[top].copy()

    other = abundance.drop(index=top, errors="ignore").sum(axis=0)
    if not other.empty:
        plot_data.loc["Other"] = other

    # Figure 1: stacked relative-abundance bar plot.
    ax = plot_data.T.plot(
        kind="bar", stacked=True, figsize=(10, 6), width=0.8
    )
    ax.set_xlabel("Sample")
    ax.set_ylabel("Relative abundance (%)")
    ax.set_title("Genus-level taxonomic composition")
    ax.legend(
        title="Genus",
        bbox_to_anchor=(1.02, 1),
        loc="upper left",
        fontsize=8
    )
    plt.xticks(rotation=45, ha="right")
    plt.tight_layout()
    plt.savefig(
        outdir / "taxonomy_barplot.png",
        dpi=300, bbox_inches="tight"
    )
    plt.close()

    # Figure 2: heatmap of the most abundant genera.
    heat = abundance.loc[top]
    fig, ax = plt.subplots(figsize=(8, max(4, len(top) * 0.3)))
    image = ax.imshow(heat.values, aspect="auto", interpolation="nearest")
    ax.set_yticks(range(len(heat.index)))
    ax.set_yticklabels(heat.index, fontsize=8)
    ax.set_xticks(range(len(sample_cols)))
    ax.set_xticklabels(sample_cols, rotation=45, ha="right")
    ax.set_xlabel("Sample")
    ax.set_ylabel("Genus")
    ax.set_title("Top genera: relative abundance (%)")
    fig.colorbar(image, ax=ax, label="Relative abundance (%)")
    fig.tight_layout()
    fig.savefig(
        outdir / "taxonomy_heatmap.png",
        dpi=300, bbox_inches="tight"
    )
    plt.close(fig)

    print(f"Samples: {len(sample_cols)}")
    print(f"Genus-level taxa: {len(abundance)}")
    print(f"Top genera plotted: {len(top)}")
    print(f"Results written to {outdir}")


if __name__ == "__main__":
    main()

