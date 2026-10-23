process SUBSAMPLE {

    tag "$sample"

    publishDir "${projectDir}/results/subsampled", mode: 'copy'

    input:
    tuple val(sample), path(read1), path(read2)

    output:
    tuple val(sample),
          path("${sample}_subsampled_1.fastq.gz"),
          path("${sample}_subsampled_2.fastq.gz"),
          emit: subsampled_reads
    
    script:
    """
    set -o pipefail

    seqtk sample -s42 ${read1} ${params.subsample_reads} | gzip -c > ${sample}_subsampled_1.fastq.gz

    seqtk sample -s42 ${read2} ${params.subsample_reads} | gzip -c > ${sample}_subsampled_2.fastq.gz
    """

    stub:
    """
    touch ${sample}_subsampled_1.fastq.gz
    touch ${sample}_subsampled_2.fastq.gz
    """
}