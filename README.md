# Gut Microbiome Shotgun Metagenomics

Reproducible **shotgun metagenomics workflow** developed in **Nextflow DSL2** to investigate gut microbiome changes associated with a legume-enriched diet in a randomized controlled trial of adults with prediabetes.

The project is being developed incrementally, with each analysis step incorporated into the workflow and documented here. The main focus is on reproducible workflow development, biological reasoning and transparent documentation of analytical decisions.


## Workflow architecture

The project is organized into two stages:

- **Part I — Quality control and taxonomic profiling (`main_I.nf`):** processes paired-end FASTQ files through raw-read quality control with FastQC and MultiQC, adapter trimming with Cutadapt, post-trimming quality control, read subsampling, taxonomic profiling with MetaPhlAn, profile merging, and taxonomic summary visualization.

- **Part II — Functional profiling (`main_II.nf`, in development):** will extend the workflow with HUMAnN to characterize microbial gene families and metabolic pathways. In progress.

---

## 1. Study background

This project uses publicly available data from the randomized controlled trial:

**A legume-enriched diet improves metabolic health in prediabetes mediated through gut microbiome: a randomized controlled trial**

* **NCBI BioProject:** `PRJNA1201249`
* **Study population:** adults with prediabetes
* **Intervention:** legume-enriched diet
* **Sample type:** fecal microbiome
* **Sequencing:** shotgun metagenomic sequencing
* **Platform:** Illumina NovaSeq 6000
* **Study duration:** 16 weeks
* **Sampling:** W0, W2, W4, W8, W12 and W16

The published study investigated changes in the gut microbiome together with metabolic and other biological measurements.

For this portfolio project, a smaller **paired longitudinal subset** was selected to develop and demonstrate a reproducible shotgun metagenomics workflow. The selected samples represent two participants with measurements at baseline (W0) and week 16 (W16).

---

## 2. Biological question

> **How does the gut microbiome change between baseline and week 16 in participants undergoing different dietary interventions in this randomized controlled trial?**

The analysis focuses on paired samples from the same participants at W0 and W16. This provides a simple longitudinal framework for exploring changes in microbial composition and function while keeping the participant as the reference unit.

The project is intended as a **portfolio-scale reproducible analysis**, rather than a complete reanalysis of the original clinical trial.

---

## 3. Samples

| Sample     | Participant | Group  | Timepoint | SRA run     |
| ---------- | ----------- | ------ | --------- | ----------- |
| FBI029_W0  | FBI029      | Group1 | W0        | SRR31800610 |
| FBI029_W16 | FBI029      | Group1 | W16       | SRR31800159 |
| FBI016_W0  | FBI016      | Group2 | W0        | SRR31800340 |
| FBI016_W16 | FBI016      | Group2 | W16       | SRR31800133 |

Metadata:

```text
samplesheet/samples.tsv
```

The complete sample metadata used during project preparation are stored in:

```text
samplesheet/shotgun_metadata.tsv
```

---

## 4. Analysis workflow

The workflow is being developed incrementally. The current planned analysis is:

```text
Raw paired-end FASTQ
        │
        ▼
      FastQC
        │
        ▼
      MultiQC
        │
        ▼
     Cutadapt
        │
        ▼
 Trimmed FASTQ
        │
        ├── FastQC
        │
        └── MultiQC
                │
                ▼
            MetaPhlAn
                │
                ▼
        Taxonomic profiling
                │
                ▼
             HUMAnN
                │
                ▼
       Functional profiling
                │
                ▼
       Downstream analysis
```

Each stage is implemented as a reproducible Nextflow component and added to the workflow as the analysis progresses.

Step-by-step analysis

Step 1 — Raw read quality control

Input: Raw paired-end shotgun metagenomic FASTQ files

The raw sequencing reads were assessed before preprocessing to evaluate sequencing quality and identify potential issues.

Tool: FastQC

FastQC was run automatically through Nextflow for all four selected biological samples, generating reports for both R1 and R2.

Output: FastQC reports

Files: results/fastqc/

Step 2 — Quality control - MultiQC

The individual FastQC reports were aggregated with MultiQC to provide a global overview of raw-read quality across the dataset.

MultiQC was implemented as a reusable Nextflow module and uses the existing FastQC .zip reports as input.

Output: results/multiqc/multiqc_report.html

## 3. Adapter identification

Before adapter trimming, the adapter sequences were identified from the sequencing data rather than assumed from the sequencing platform alone.

The study used paired-end Illumina NovaSeq 6000 sequencing. FastQC detected low-level adapter contamination in the raw reads. To investigate the adapter sequences, the expected Illumina adapter sequences were searched directly in the FASTQ files.

The following sequences were detected in the selected sample (`FBI029_W0`):

| Read | Adapter sequence                    | Matches |
| ---- | ----------------------------------- | ------: |
| R1   | `AGATCGGAAGAGCACACGTCTGAACTCCAGTCA` |   7,477 |
| R2   | `AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT` |   7,532 |

The presence of these sequences in both reads supports their use as the adapter sequences for the paired-end trimming step with Cutadapt.

This provides a data-supported basis for the adapter configuration used in the next workflow step.

### 4. Adapter trimming

The identified adapter sequences were removed from the paired-end reads using Cutadapt.
Cutadapt was implemented as a reusable Nextflow module and processes all samples automatically.
A JSON report was generated for each sample to retain the trimming statistics.

Input: file_1.fastq.gz file_2.fastq.gz
Output: results/trimmed/

A summary table was generated automatically from the Cutadapt reports:
Output: results/trimmed/cutadapt_summary.tsv

The trimming step was completed successfully for all four paired-end samples. Adapter sequences were detected in approximately 2.4–2.8% of reads, while the majority of reads did not contain detectable adapter sequence.

| Sample     | Input reads | R1 adapter reads (%) | R2 adapter reads (%) |
| ---------- | ----------: | -------------------: | -------------------: |
| FBI029_W0  |  24,354,291 |      596,811 (2.45%) |      583,191 (2.39%) |
| FBI029_W16 |  23,684,666 |      667,974 (2.82%) |      645,015 (2.72%) |
| FBI016_W0  |  24,815,969 |      721,920 (2.91%) |      697,592 (2.81%) |
| FBI016_W16 |  25,160,519 |      884,940 (3.52%) |      878,403 (3.49%) |


## 5. Quality control after trimming

FastQC was run again on the trimmed reads to assess read quality after adapter removal.
The FastQC results were aggregated with MultiQC to provide a global overview of the processed reads.

Both steps were implemented as reusable Nextflow modules.

Output: results/fastqc_trimmed/

Output: results/multiqc_trimmed/multiqc_report.html

The post-trimming QC was used to verify that the processed reads were suitable for downstream taxonomic profiling.

## 6. Read subsampling for taxonomic profiling

A subset of the trimmed reads was generated before taxonomic profiling to enable an efficient local benchmark with MetaPhlAn.

Because the complete dataset contains approximately 24 million paired-end reads per sample, a full MetaPhlAn mapping run was computationally demanding on the local system. A subset of 1 million paired reads per sample was therefore selected for this portfolio-scale analysis.

Subsampling was implemented as a reusable Nextflow module using `seqtk`. The same random seed (`42`) was used for R1 and R2 to maintain paired-end correspondence.

Input: `results/trimmed/<sample>_trimmed_1.fastq.gz`  
       `results/trimmed/<sample>_trimmed_2.fastq.gz`

Output: `results/subsampled/<sample>_subsampled_1.fastq.gz`  
        `results/subsampled/<sample>_subsampled_2.fastq.gz`

The subsampling process was applied automatically to all four samples.

The 1 million-read subset was used as a computational benchmark for downstream taxonomic profiling rather than as a replacement for the complete sequencing depth.

## 7. Taxonomic profiling and profile merging with MetaPhlAn

MetaPhlAn was used to estimate the taxonomic composition of the shotgun metagenomic samples from the subsampled paired-end reads.

The workflow supports two execution modes, controlled by the `use_existing_merge` parameter in `nextflow.config`:

* When `use_existing_merge = false`, the `METAPHLAN` module runs taxonomic profiling on the subsampled reads and produces an individual profile for each sample.
* When `use_existing_merge = true`, the workflow reuses existing MetaPhlAn profiles from `results/metaphlan/` and skips the profiling step.

In both modes, the individual profiles are passed to the `METAPHLAN_MERGE` module, defined in `modules/metaphlan_merge/main_metaphlan_merge.nf`. This module combines the profiles using `merge_metaphlan_tables.py` and generates a single table containing taxa as rows and samples as columns.

This design separates taxonomic profiling from profile merging, allowing the downstream steps to be rerun without repeating the computationally intensive profiling process.

**Input:** `results/metaphlan/<sample>_profile.txt` when reusing existing profiles.

**Output:** `results/metaphlan_merge/metaphlan_merged.tsv`

## 8. Genus-level taxonomic summary and visualization

The merged MetaPhlAn table is passed directly to the `TAXONOMY_SUMMARY` module, defined in `modules/taxonomy_summary/main_taxonomy_summary.nf`.

This module executes the Python script `scripts/taxonomy_summary.py`, which extracts genus-level entries from the merged table and summarizes their relative abundances across samples. Restricting the summary to genus-level entries avoids mixing different taxonomic ranks in the same visualization.

The analysis generates three outputs:

* `results/taxonomy_summary/taxonomy_summary.tsv`: genus-level relative abundance table.
* `results/taxonomy_summary/taxonomy_barplot.png`: stacked bar plot showing the 15 genera with the highest mean relative abundance across samples. The remaining genera are grouped as `Other`.
* `results/taxonomy_summary/taxonomy_heatmap.png`: heatmap comparing the relative abundances of the same 15 genera across samples.

### Genus-level taxonomic composition

The stacked bar plot provides an overview of the relative contribution of the most abundant genera in each sample.

![Genus-level taxonomic composition](results/taxonomy_summary/taxonomy_barplot.png)

### Comparison of the most abundant genera

The heatmap complements the bar plot by highlighting differences in the relative abundance patterns of the selected genera across samples.

![Heatmap of the top 15 genera](results/taxonomy_summary/taxonomy_heatmap.png)

These visualizations provide an exploratory overview of microbial community composition and can help identify taxa for further investigation.

The current analysis includes four samples: `FBI029_W0`, `FBI029_W16`, `FBI016_W0`, and `FBI016_W16`. Comparisons between baseline (W0) and week 16 (W16) are descriptive and exploratory. This small subset is insufficient to establish robust diet-associated effects or statistical significance.

The modular design allows the taxonomic summary and visualizations to be regenerated from the merged table without rerunning the upstream MetaPhlAn profiling step.













## Reference databases

### MetaPhlAn

- **Version:** 4.1.1
- **Database index:** `mpa_vJun23_CHOCOPhlAnSGB_202403`
- **Purpose:** Taxonomic profiling of shotgun metagenomic reads.
- **Storage:** External storage, outside the Git repository.
- **Configuration:** The database directory and index are specified in the local Nextflow configuration.

### HUMAnN

| Database | Purpose | Storage |
|---|---|---|
| ChocoPhlAn, full | Nucleotide-level search against microbial pangenomes | External storage |
| UniRef50, full annotated | Protein-level search and functional annotation | External storage |

Example local configuration paths:

```groovy
params {
    metaphlan_db = '/path/to/metaphlan_database_humann'
    metaphlan_index = 'mpa_vJun23_CHOCOPhlAnSGB_202403'

    humann_chocophlan_db = '/path/to/humann_database/chocophlan'
    humann_uniref_db = '/path/to/humann_database/uniref'
}
```

Replace these example paths with the actual database locations on your machine before running the workflow.

Storage and reproducibility

Keep reference databases outside the Git repository.

Do not commit FASTQ files, database files, or large intermediate outputs.

Store machine-specific paths in a local Nextflow configuration rather than hard-coding them into reusable workflow modules.

Ensure sufficient free disk space for reference databases, input FASTQs, temporary files, and workflow outputs.

Use persistent storage for databases and final outputs, and appropriate temporary storage for intermediate files.

Record software versions and database index names for each analysis to support reproducibility.

Reproducibility note: The current configuration uses MetaPhlAn 4.1.1 with mpa_vJun23_CHOCOPhlAnSGB_202403 and HUMAnN 3.9 with the full ChocoPhlAn and UniRef50 databases. Changing software versions or database releases may change taxonomic and functional profiling results. Database paths are specific to the local environment and must be adapted when running the workflow on another machine.





