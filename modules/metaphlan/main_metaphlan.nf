process METAPHLAN {

    maxForks 1

    tag "$sample"

    publishDir "${projectDir}/results/metaphlan", mode: 'copy'

    input:
    tuple val(sample), path(read1), path(read2)

    output:
    tuple val(sample), path("${sample}_profile.txt"), emit: profile
    path "${sample}_mapout.bz2", emit: mapout

    script:
    """
    metaphlan \
        ${read1},${read2} \
        --input_type fastq \
        --db_dir ${params.metaphlan_db} \
        --index ${params.metaphlan_index} \
        --offline \
        --nproc ${params.metaphlan_threads} \
        --mapout ${sample}_mapout.bz2 \
        -o ${sample}_profile.txt
    """
        stub:
    """
    touch ${sample}_profile.txt
    touch ${sample}_mapout.bz2
    """
}
