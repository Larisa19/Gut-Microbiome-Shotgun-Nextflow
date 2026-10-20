nextflow.enable.dsl=2

include { MULTIQC } from './modules/multiqc/main_multiqc.nf'
include { CUTADAPT; SUMMARY_CUTADAPT } from './modules/cutadapt/main_cutadapt'
include { FASTQC_TRIMMED } from './modules/fastqc_trimmed/main_fastqc_trimmed.nf'
include { MULTIQC_TRIMMED } from './modules/multiqc_trimmed/main_multiqc_trimmed.nf'
include { SUBSAMPLE } from './modules/subsample/main_subsample.nf'
include { METAPHLAN } from './modules/metaphlan/main_metaphlan.nf'
include { METAPHLAN_MERGE } from './modules/metaphlan_merge/main_metaphlan_merge.nf'
include { TAXONOMY_SUMMARY } from './modules/taxonomy_summary/main_taxonomy_summary.nf'

process FASTQC_RAW {

    tag "$sample"

    publishDir "${projectDir}/results/fastqc", mode: 'copy'

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


workflow {

    // 1. Read raw paired-end FASTQ files
    samples = channel
        .fromPath('data/raw/*_1.fastq.gz')
        .map { read1 ->

            def sample = read1.name.replaceFirst(/_1\.fastq\.gz$/, '')
            def read2 = file("data/raw/${sample}_2.fastq.gz")

            tuple(sample, read1, read2)
        }
    // 2. Raw read quality control
    FASTQC_RAW(samples)

        MULTIQC(
        FASTQC_RAW.out.fastqc_zip
            .map { sample, zip -> zip }
            .collect()
    )
    // 3. Adapter trimming
    samples_cutadapt = channel
        .fromPath('data/raw/*_1.fastq.gz')
        .map { read1 ->
            def sample = read1.name.replaceFirst(/_1\.fastq\.gz$/, '')
            def read2 = file("data/raw/${sample}_2.fastq.gz")

            tuple(sample, read1, read2)
        }

    CUTADAPT(samples_cutadapt)

    SUMMARY_CUTADAPT(
        CUTADAPT.out.cutadapt_report.collect()
    )
    
    // 4. Quality control after trimming
    FASTQC_TRIMMED(
        CUTADAPT.out.trimmed_reads
    )

    trimmed_fastqc = FASTQC_TRIMMED.out.fastqc_zip
    .map { sample, zip -> zip }
    .collect()

    MULTIQC_TRIMMED(trimmed_fastqc)
    
    // 5. Subsample reads
    SUBSAMPLE(
        CUTADAPT.out.trimmed_reads
    )

    
    // 6. MetaPhlAn: reuse existing profiles or run MetaPhlAn
    if (params.use_existing_merge) {

        profiles = channel.fromPath(
            'results/metaphlan/*_profile.txt',
            checkIfExists: true
        )

    } else {

        METAPHLAN(
            SUBSAMPLE.out.subsampled_reads
        )

        profiles = METAPHLAN.out.profile
            .map { sample, profile -> profile }
    }

    // 7. Merge MetaPhlAn profiles
    
    METAPHLAN_MERGE(
        profiles.collect()
    )

    merged_table = METAPHLAN_MERGE.out.merged_table

    // 8. Taxonomic summary
    TAXONOMY_SUMMARY(merged_table)
}