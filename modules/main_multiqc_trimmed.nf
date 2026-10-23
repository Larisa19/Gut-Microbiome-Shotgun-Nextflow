process MULTIQC_TRIMMED {

    tag "Trimmed read QC"

    publishDir "${projectDir}/results/multiqc_trimmed", mode: 'copy'

    input:
    path fastqc_results

    output:
    path "multiqc_report.html"

    script:
    """
    multiqc . \
        --force \
        --filename multiqc_report.html
    """
}
