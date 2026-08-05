#!/usr/bin/env python3
"""
Merge individual CNV depth matrices into a single sample x region matrix.

This script:
1. Reads individual cnv_per_region.tsv files for each sample
2. Extracts region and log2cnv columns
3. Pivots to create sample x region matrix (transposed)
4. Formats columns as "cnv_depth__{region_name}"
5. Outputs merged CNV depth matrix

Input: One or more cnv_per_region.tsv files (tab-separated values)
Output: cnv_depth_matrix.tsv (sample x region matrix)
"""

import argparse
import pandas as pd


def merge_cnv_matrices(input_files, sample_names, output_file):
    """Merge individual CNV depth matrices into a single matrix."""
    dfs = []

    for sample, tsv in zip(sample_names, input_files):
        df = pd.read_csv(tsv, sep="\t")

        # Original: log2cnv depth ratios
        df_log2 = df[["region", "log2cnv"]].copy()
        df_log2.columns = ["region", sample]
        df_log2 = df_log2.set_index("region")
        dfs.append(df_log2)

    # Log2 matrix (original format)
    matrix = pd.concat(dfs, axis=1).T
    matrix.index.name = "sample"
    matrix.columns = [f"cnv_depth__{c}" for c in matrix.columns]
    matrix.to_csv(output_file, sep="\t")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Merge individual CNV depth matrices into a single sample x region matrix"
    )
    parser.add_argument(
        "--input-files",
        nargs="+",
        required=True,
        help="List of cnv_per_region.tsv files to merge"
    )
    parser.add_argument(
        "--sample-names",
        nargs="+",
        required=True,
        help="List of sample names corresponding to input files"
    )
    parser.add_argument(
        "--output",
        required=True,
        help="Output merged CNV depth matrix file"
    )

    args = parser.parse_args()

    if len(args.input_files) != len(args.sample_names):
        raise ValueError(
            f"Number of input files ({len(args.input_files)}) "
            f"does not match number of sample names ({len(args.sample_names)})"
        )

    merge_cnv_matrices(args.input_files, args.sample_names, args.output)