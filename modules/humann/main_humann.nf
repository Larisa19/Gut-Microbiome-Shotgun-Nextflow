process HUMANN {

    tag "${sample}"

    cpus 2
    memory '12 GB'
    time '48h'

    publishDir "${projectDir}/results/humann", mode: 'copy'

    input:
    tuple val(sample),
          path(read1),
          path(read2),
          path(profile)

    output:
    tuple val(sample),
          path("${sample}_humann"),
          emit: humann_results

    script:
    """
    set -euo pipefail

    # Concatenate paired-end reads; this does not merge read overlaps.
    gzip -cd "${read1}" "${read2}" | gzip -c > "${sample}.fastq.gz"

    humann \
        --input "${sample}.fastq.gz" \
        --input-format fastq.gz \
        --output "${sample}_humann" \
        --threads ${task.cpus} \
        --nucleotide-database "${params.humann_chocophlan_db}" \
        --protein-database "${params.humann_uniref_db}" \
        --taxonomic-profile "${profile}" \
        --remove-temp-output
    """
}