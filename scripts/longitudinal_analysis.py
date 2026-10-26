
#!/usr/bin/env python3
"""
Longitudinal analysis of a shotgun metagenomics project.

Compares W0 vs W16 for participants FBI016 and FBI029.
- Normalizes HUMAnN pathway abundance to CPM using humann_renorm_table.
- Calculates paired changes in genus abundance and pathway abundance.
- Generates TSV tables and paired line plots.
- Performs descriptive analysis only; no statistical significance tests.

Run from any directory:
    python scripts/longitudinal_analysis.py
"""

from pathlib import Path
import shutil
import subprocess
import sys

import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


# ---------------------------------------------------------
# 1. Paths and sample definitions
# ---------------------------------------------------------

PROJECT_ROOT = Path(__file__).resolve().parent.parent

TAXONOMY_FILE = (
    PROJECT_ROOT / "results" / "taxonomy_summary" / "taxonomy_summary.tsv"
)

HUMANN_FILE = (
    PROJECT_ROOT / "results" / "humann" / "pathabundance_merged.tsv"
)

OUTPUT_DIR = PROJECT_ROOT / "results" / "longitudinal"
NORMALIZED_FILE = (
    PROJECT_ROOT / "results" / "longitudinal" / "pathabundance_cpm.tsv"
)

PARTICIPANTS = {
    "FBI016": ("FBI016_W0", "FBI016_W16"),
    "FBI029": ("FBI029_W0", "FBI029_W16"),
}

TOP_N = 15


# ---------------------------------------------------------
# 2. Helper functions
# ---------------------------------------------------------

def require_file(path):
    """Stop with a clear message if an input file is missing."""
    if not path.is_file():
        raise FileNotFoundError(
            f"Required file not found:\n{path}\n"
            "Check that the corresponding workflow has finished "
            "and that the file is in the expected location."
        )


def find_humann_renorm():
    """Find the HUMAnN normalization executable."""
    executable = shutil.which("humann_renorm_table")

    if executable:
        return executable

    # Common location for the HUMAnN Conda environment on this Mac.
    fallback = Path("/opt/miniconda3/envs/humann/bin/humann_renorm_table")

    if fallback.is_file():
        return str(fallback)

    raise FileNotFoundError(
        "Could not find humann_renorm_table.\n"
        "Activate the HUMAnN environment or install/configure HUMAnN "
        "so its executable is available."
    )


def read_table(path):
    """Read a tab-separated table and clean column names."""
    df = pd.read_csv(path, sep="\t")

    if df.empty:
        raise ValueError(f"The input table is empty: {path}")

    df.columns = [str(c).strip().lstrip("#").strip() for c in df.columns]

    if df.columns.duplicated().any():
        raise ValueError(f"Duplicate column names found in {path}")

    return df


def validate_samples(df, required_samples, table_name):
    """Verify that every expected sample is present."""
    missing = [s for s in required_samples if s not in df.columns]

    if missing:
        raise ValueError(
            f"{table_name} is missing expected sample columns: {missing}\n"
            f"Available columns: {list(df.columns)}"
        )


def make_numeric(df, columns, table_name):
    """Convert abundance columns to numeric values."""
    for col in columns:
        df[col] = pd.to_numeric(df[col], errors="coerce")

    if df[columns].isna().any().any():
        raise ValueError(
            f"Missing or non-numeric abundance values found in {table_name}."
        )

    if (df[columns] < 0).any().any():
        raise ValueError(f"Negative abundance values found in {table_name}.")


def paired_changes(df, feature_col, sample_columns):
    """
    Create one row per feature and participant, with W0, W16,
    absolute difference and log2 fold change when defined.
    """
    records = []

    for _, row in df.iterrows():
        feature = row[feature_col]

        for participant, (w0_sample, w16_sample) in PARTICIPANTS.items():
            w0 = float(row[w0_sample])
            w16 = float(row[w16_sample])

            if w0 > 0 and w16 > 0:
                log2_fc = __import__("math").log2(w16 / w0)
            else:
                # Undefined when either value is zero.
                log2_fc = float("nan")

            records.append({
                "Feature": feature,
                "Participant": participant,
                "W0": w0,
                "W16": w16,
                "Difference_W16_minus_W0": w16 - w0,
                "Log2_fold_change": log2_fc,
            })

    return pd.DataFrame(records)


def plot_paired_features(
    df,
    feature_col,
    sample_columns,
    ranking_column,
    title,
    ylabel,
    output_path,
):
    """Plot W0-to-W16 changes for the top features in each table."""
    means = df.set_index(feature_col)[sample_columns].mean(axis=1)
    top_features = means.nlargest(TOP_N).index.tolist()

    if not top_features:
        print(f"Skipping empty plot: {output_path.name}")
        return

    fig, axes = plt.subplots(
        nrows=len(top_features),
        ncols=1,
        figsize=(9, max(3, 2.0 * len(top_features))),
        squeeze=False,
    )

    for ax, feature in zip(axes.flatten(), top_features):
        feature_row = df.loc[df[feature_col] == feature].iloc[0]

        for participant, (w0_sample, w16_sample) in PARTICIPANTS.items():
            values = [
                float(feature_row[w0_sample]),
                float(feature_row[w16_sample]),
            ]
            ax.plot(
                [0, 1],
                values,
                marker="o",
                linewidth=1.5,
                label=participant,
            )

        ax.set_title(str(feature), loc="left", fontsize=9)
        ax.set_xticks([0, 1], ["W0", "W16"])
        ax.set_ylabel(ylabel, fontsize=8)
        ax.grid(axis="y", alpha=0.3)

    axes[0, 0].legend(loc="best", fontsize=8)
    fig.suptitle(title, y=1.002)
    fig.tight_layout()
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)

    print(f"Created: {output_path}")


# ---------------------------------------------------------
# 3. Check input files and create output directory
# ---------------------------------------------------------

def main():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    require_file(TAXONOMY_FILE)
    require_file(HUMANN_FILE)

    all_samples = [
        sample
        for pair in PARTICIPANTS.values()
        for sample in pair
    ]

    print("Project:", PROJECT_ROOT)
    print("Output directory:", OUTPUT_DIR)


    # -----------------------------------------------------
    # 4. Taxonomic longitudinal analysis
    # -----------------------------------------------------

    print("\n[1/4] Analysing genus abundances...")

    taxonomy = read_table(TAXONOMY_FILE)

    if "Genus" not in taxonomy.columns:
        raise ValueError(
            f"Expected a 'Genus' column in {TAXONOMY_FILE}"
        )

    validate_samples(taxonomy, all_samples, "Taxonomy table")
    make_numeric(taxonomy, all_samples, "Taxonomy table")

    if taxonomy["Genus"].duplicated().any():
        raise ValueError("Duplicate genus names found in taxonomy table.")

    taxonomy = taxonomy.dropna(subset=["Genus"])
    taxonomy = taxonomy[taxonomy["Genus"].astype(str).str.strip() != ""]

    taxonomy_changes = paired_changes(
        taxonomy, "Genus", all_samples
    )

    taxonomy_changes.to_csv(
        OUTPUT_DIR / "taxonomy_changes.tsv",
        sep="\t",
        index=False,
    )

    taxonomy.to_csv(
        OUTPUT_DIR / "taxonomy_abundances.tsv",
        sep="\t",
        index=False,
    )

    plot_paired_features(
        taxonomy,
        feature_col="Genus",
        sample_columns=all_samples,
        ranking_column=all_samples,
        title="Longitudinal genus abundance: W0 vs W16",
        ylabel="Abundance (%)",
        output_path=OUTPUT_DIR / "taxonomy_paired_changes.png",
    )


    # -----------------------------------------------------
    # 5. Normalize HUMAnN pathway abundances to CPM
    # -----------------------------------------------------

    print("\n[2/4] Normalizing HUMAnN pathway abundances to CPM...")

    renorm_executable = find_humann_renorm()

    subprocess.run(
        [
            renorm_executable,
            "--input", str(HUMANN_FILE),
            "--output", str(NORMALIZED_FILE),
            "--units", "cpm",
            "--mode", "community",
        ],
        check=True,
    )

    require_file(NORMALIZED_FILE)
    print(f"Created normalized table: {NORMALIZED_FILE}")


    # -----------------------------------------------------
    # 6. Functional longitudinal analysis
    # -----------------------------------------------------

    print("\n[3/4] Analysing pathway abundances...")

    pathways = read_table(NORMALIZED_FILE)

    if "Pathway" not in pathways.columns:
        raise ValueError(
            f"Expected a 'Pathway' column in {NORMALIZED_FILE}"
        )

    validate_samples(pathways, all_samples, "HUMAnN pathway table")
    make_numeric(pathways, all_samples, "HUMAnN pathway table")

    # Keep community-level pathways only:
    # exclude unmapped/unintegrated rows and stratified rows
    # whose names contain a vertical bar.
    names = pathways["Pathway"].astype(str)

    excluded = (
        names.isin(["UNMAPPED", "UNINTEGRATED"])
        | names.str.startswith("UNINTEGRATED|")
        | names.str.contains("|", regex=False)
    )

    pathways = pathways.loc[~excluded].copy()

    if pathways.empty:
        raise ValueError(
            "No community-level pathways remain after filtering."
        )

    if pathways["Pathway"].duplicated().any():
        raise ValueError("Duplicate pathway names found in HUMAnN table.")

    pathway_changes = paired_changes(
        pathways, "Pathway", all_samples
    )

    pathway_changes.to_csv(
        OUTPUT_DIR / "pathway_changes.tsv",
        sep="\t",
        index=False,
    )

    pathways.to_csv(
        OUTPUT_DIR / "pathway_abundances_cpm.tsv",
        sep="\t",
        index=False,
    )

    plot_paired_features(
        pathways,
        feature_col="Pathway",
        sample_columns=all_samples,
        ranking_column=all_samples,
        title="Longitudinal pathway abundance: W0 vs W16",
        ylabel="Abundance (CPM)",
        output_path=OUTPUT_DIR / "pathway_paired_changes.png",
    )


    # -----------------------------------------------------
    # 7. Summary and completion message
    # -----------------------------------------------------

    print("\n[4/4] Analysis complete.")
    print("\nOutput files:")

    for path in sorted(OUTPUT_DIR.iterdir()):
        if path.is_file():
            print(f" - {path.relative_to(PROJECT_ROOT)}")

    print(
        "\nInterpretation note: these are descriptive paired comparisons "
        "of two participants. They do not establish a statistically "
        "significant or causal diet effect."
    )


if __name__ == "__main__":
    try:
        main()
    except (FileNotFoundError, ValueError, subprocess.CalledProcessError) as exc:
        print(f"\nERROR: {exc}", file=sys.stderr)
        sys.exit(1)


