process FASTQC_TRIMMED {

    tag "$sample"

    publishDir "${projectDir}/results/fastqc_trimmed", mode: 'copy'

    input:
    tuple val(sample), path(read1), path(read2)

    output:
    tuple val(sample), path("*_fastqc.zip"), emit: fastqc_zip
    path "*_fastqc.html"

    script:
    """
    fastqc ${read1} ${read2}
    """
}