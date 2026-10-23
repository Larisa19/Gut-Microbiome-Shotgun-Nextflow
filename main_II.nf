nextflow.enable.dsl = 2

include { HUMANN } from './modules/humann/main_humann.nf'
include { HUMANN_MERGE } from './modules/humann/main_humann_merge.nf'

workflow {

    reads_ch = Channel.fromFilePairs(
        "${projectDir}/results/subsampled/*_{1,2}.fastq.gz",
        checkIfExists: true
    )

    humann_input_ch = reads_ch.map { sample, reads ->
        def sample_id = sample.replaceFirst(/_subsampled$/, '')

        def profile = file(
            "${projectDir}/results/metaphlan/${sample_id}_profile.txt",
            checkIfExists: true
        )

                println "INPUT: ${sample_id} | R1=${reads[0].name} | R2=${reads[1].name} | PROFILE=${profile.name}"

        tuple(sample_id, reads[0], reads[1], profile)
    }

    HUMANN(humann_input_ch)
    HUMANN_MERGE(HUMANN.out.humann_results.map { sample, humann_dir -> humann_dir }.collect())    
}