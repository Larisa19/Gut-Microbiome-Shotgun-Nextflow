
process TAXONOMY_SUMMARY {

    tag "Genus-level taxonomic summary"

    publishDir "${projectDir}/results/taxonomy_summary", mode: 'copy'

    input:
    path merged_table

    output:
    path "taxonomy_summary.tsv", emit: summary
    path "taxonomy_barplot.png", emit: barplot
    path "taxonomy_heatmap.png", emit: heatmap

    script:
    """
    python "${projectDir}/scripts/taxonomy_summary.py" \
        --input "${merged_table}" \
        --outdir "." \
        --top-n 15
    """
}

