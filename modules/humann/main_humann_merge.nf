
process HUMANN_MERGE {

    tag 'all_samples'

    cpus 1

    publishDir "${projectDir}/results/humann", mode: 'copy'

    input:
    path humann_dirs

    output:
    path 'genefamilies_merged.tsv', emit: genefamilies
    path 'pathabundance_merged.tsv', emit: pathabundance

    script:
    """
    set -euo pipefail

    mkdir genefamilies pathabundance

    for sample_dir in *_humann; do
        cp "\${sample_dir}"/*_genefamilies.tsv genefamilies/
        cp "\${sample_dir}"/*_pathabundance.tsv pathabundance/
    done

    humann_join_tables \
        --input genefamilies \
        --output genefamilies_merged.tsv \
        --file_name genefamilies

    humann_join_tables \
        --input pathabundance \
        --output pathabundance_merged.tsv \
        --file_name pathabundance
    """
}
