process CUTADAPT {

    maxForks 1
    
    tag "$sample"

    publishDir "${projectDir}/results/trimmed", mode: 'copy'

    input:
    tuple val(sample), path(read1), path(read2)

    output:
    tuple val(sample),
          path("${sample}_trimmed_1.fastq.gz"),
          path("${sample}_trimmed_2.fastq.gz"),
          emit: trimmed_reads
    path "${sample}_cutadapt.json", emit: cutadapt_report

    script:
    """
    cutadapt \
        -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
        -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
        -o ${sample}_trimmed_1.fastq.gz \
        -p ${sample}_trimmed_2.fastq.gz \
        --json ${sample}_cutadapt.json \
        ${read1} ${read2}
    """
}
process SUMMARY_CUTADAPT {

    publishDir "${projectDir}/results/trimmed", mode: 'copy'

    input:
    path reports

    output:
    path "cutadapt_summary.tsv"

    script:
    """
    python ${projectDir}/scripts/summarize_cutadapt.py
    """
}