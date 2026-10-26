process METAPHLAN_MERGE {

    tag "Merge MetaPhlAn profiles"

    publishDir "${projectDir}/results/metaphlan", mode: 'copy'

    input:
    path profiles

    output:
    path "metaphlan_merged.tsv", emit: merged_table

    script:
    """
    merge_metaphlan_tables.py ${profiles} > metaphlan_merged.tsv
    """
}