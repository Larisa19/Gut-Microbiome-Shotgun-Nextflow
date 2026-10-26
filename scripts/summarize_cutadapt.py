import json
import glob
import os

json_files = sorted(glob.glob("*_cutadapt.json"))

with open("cutadapt_summary.tsv", "w") as out:
    out.write(
        "Sample\tInput_reads\tR1_adapter_reads\tR1_adapter_%\t"
        "R2_adapter_reads\tR2_adapter_%\n"
    )

    for file in json_files:
        with open(file) as f:
            data = json.load(f)

        sample = os.path.basename(file).replace("_cutadapt.json", "")
        input_reads = data["read_counts"]["input"]

        r1 = data["adapters_read1"][0]["total_matches"]
        r2 = data["adapters_read2"][0]["total_matches"]

        r1_pct = 100 * r1 / input_reads
        r2_pct = 100 * r2 / input_reads

        out.write(
            f"{sample}\t{input_reads}\t{r1}\t{r1_pct:.2f}\t"
            f"{r2}\t{r2_pct:.2f}\n"
        )
