import pandas as pd 
import yaml
import json

# Load sample names from the config file
with open("config/config.yaml") as f:
    config = yaml.safe_load(f)
samples = config["samples"]["name"]

# Load multiqc results
with open("results/00_QC/multiqc/multiqc_data/multiqc_data.json", "r") as f:
    data = json.load(f)
multiqc_df = pd.DataFrame(data["report_general_stats_data"]["fastp"]).transpose()

# Load coverage results
covs = []
for sample in samples:
    cov = pd.read_csv(f"results/00_QC/coverage/{sample}/coverage.txt", sep="\t", header=0)
    cov["sample"] = sample
    covs.append(cov)
coverage_stats_df = pd.concat(covs, ignore_index=True)
coverages = coverages = coverage_stats_df.groupby("sample")["coverage"].mean()

# make sure all samples are in the multiqc file
missing_samples = [s for s in samples if s not in multiqc_df.index]
if missing_samples:
    print("Not all samples in MultiQC file!")
    raise ValueError("Missing samples in MultiQC file")
    
# make sure all samples are in the coverages file
missing_samples_coverage = [s for s in samples if s not in coverages.index]
if missing_samples_coverage:
    print("Not all samples in coverage files!")
    raise ValueError("Missing samples in coverage files")

# Define filtering thresholds
filters = {
    "after_filtering_q30_rate": 0.9,
    "expected_gc_content": 0.67,
    "percent_filtered_reads": 0.75,
    "coverage_min": 0.01
}

# Create masks (will be false for samples we want to keep)
q30_filter = ~(multiqc_df["after_filtering_q30_rate"] > filters["after_filtering_q30_rate"])
gc_filter = ~(multiqc_df["after_filtering_gc_content"].between(filters["expected_gc_content"]-0.03, filters["expected_gc_content"] + 0.03))
percent_passed_filter = ~(multiqc_df["pct_surviving"] > filters["percent_filtered_reads"])
coverage_filter = ~(coverages > filters["coverage_min"])

passed_filter = q30_filter | gc_filter | percent_passed_filter | coverage_filter

# filter the multiqc data frame
multiqc_df = multiqc_df[~passed_filter]

# passed samples
sample_passed = multiqc_df.index.values

pd.Series(sample_passed).to_csv("results/00_QC/passed_samples.csv", index=False, header=False)

