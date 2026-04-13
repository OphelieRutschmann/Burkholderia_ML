import pandas as pd 
import yaml

# Load sample names from the config file
with open("config/config.yaml") as f:
    config = yaml.safe_load(f)
samples = config["samples"]["name"]

# Load multiqc results
multiqc_df = pd.read_csv("results/00_QC/multiqc/multiqc_data/multiqc_fastp.txt", sep="\t")

# make sure all samples are in the multiqc file
missing_samples = [s for s in samples if s not in multiqc_df["sample"].values]
if missing_samples:
    print("Not all samples in MultiQC file!")
    raise ValueError("Missing samples in MultiQC file")

# Define filtering thresholds
filters = {
    "after_filtering_q30_rate": 0.9,
    "after_filtering_gc_content": 50,
    "filtering_result_low_quality_reads": 90,
    "duplication_rate": None
}

# Create masks (will be false for samples we want to keep)
q30_filter = ~(multiqc_df["after_filtering_q30_rate"] > filters["after_filtering_q30_rate"])

gc_filter = ~(multiqc_df["after_filtering_gc_content"].between(
    filters["after_filtering_gc_content"] - 10,
    filters["after_filtering_gc_content"] + 10
))

# will be true for any row that we should remove (that failed at least one of the quality checks)
mask = (q30_filter | gc_filter)

print(mask.sum(), "samples will be removed")

# filter the multiqc data frame
multiqc_df = multiqc_df[~mask]

# passed samples
sample_passed = multiqc_df["sample"]

sample_passed.to_csv("results/00_QC/passed_samples.txt", index=False, header=False)
