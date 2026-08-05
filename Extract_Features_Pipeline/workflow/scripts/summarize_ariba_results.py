#!/usr/bin/env python3
"""
Summarize ARIBA resistance gene detection results into a CSV matrix.

Input: ARIBA report.tsv files from multiple samples
Output: CSV with samples x resistance genes matrix
"""

import pandas as pd
import os
from pathlib import Path

# Get parameters from Snakemake
min_identity = float(snakemake.params.min_identity)
min_coverage = float(snakemake.params.min_coverage)

# Parse all ARIBA reports
detected_genes = {}
all_genes = set()

for report_file in snakemake.input.reports:
    sample_name = Path(report_file).parent.name
    detected_genes[sample_name] = {}

    if not os.path.exists(report_file):
        print(f"WARNING: Report not found for {sample_name}: {report_file}")
        continue

    try:
        df = pd.read_csv(report_file, sep="\t")

        for _, row in df.iterrows():
            gene = row.get('ref_name', 'unknown')
            identity = float(row.get('query_identity', 0))
            coverage = float(row.get('ref_coverage', 0))
            status = row.get('status', 'match')

            # Filter by quality thresholds
            if identity >= min_identity and coverage >= min_coverage:
                detected_genes[sample_name][gene] = status
                all_genes.add(gene)

    except Exception as e:
        print(f"ERROR parsing {report_file}: {e}")
        continue

# Build summary matrix
sorted_genes = sorted(all_genes)
summary_data = []

for sample_name in detected_genes.keys():
    row = {'sample': sample_name}
    sample_genes = detected_genes[sample_name]

    for gene in sorted_genes:
        if gene in sample_genes:
            row[gene] = sample_genes[gene]
        else:
            row[gene] = "absent"

    summary_data.append(row)

# Write summary CSV
summary_df = pd.DataFrame(summary_data)
summary_df.to_csv(snakemake.output.summary, index=False)

print(f"Summary: {len(summary_df)} samples, {len(sorted_genes)} unique genes")
