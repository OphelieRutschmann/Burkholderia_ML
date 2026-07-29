#!/usr/bin/env python3
"""
Normalize CNV depth values to log2 ratios per chromosome.

This script:
1. Reads mosdepth region depth output (gzipped BED format)
2. Calculates per-chromosome median depth using only CDS and gene features
3. Computes CNV ratio (depth / chromosome_median)
4. Calculates log2 CNV values (log2(cnv_ratio))
5. Outputs normalized CNV matrix per region

Input: mosdepth.regions.bed.gz (gzip-compressed BED file)
Output: cnv_per_region.tsv (tab-separated values)
"""

import argparse
import numpy as np
import pandas as pd


def normalize_cnv_depth(input_file, output_file):
    """Normalize CNV depth to log2 ratios per chromosome."""
    df = pd.read_csv(
        input_file,
        sep="\t",
        compression="gzip",
        header=None,
        names=["chrom", "start", "end", "region", "depth"]
    )

    # Use only CDS and gene features to estimate a stable per-chromosome coverage
    # Don't consider promoter regions
    mask = (
        df["region"].str.contains("__CDS", na=False)
        | df["region"].str.contains("__gene", na=False)
    )
    medians = df.loc[mask].groupby("chrom")["depth"].median().to_dict()

    df["chrom_median"] = df["chrom"].map(medians)
    df["cnv_ratio"] = df["depth"] / df["chrom_median"]
    # Avoid -inf when depth == 0
    df["log2cnv"] = np.log2(df["cnv_ratio"].replace(0, np.nan))

    df[["region", "chrom", "depth", "chrom_median", "cnv_ratio", "log2cnv"]].to_csv(
        output_file,
        sep="\t",
        index=False
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Normalize CNV depth to log2 ratios per chromosome"
    )
    parser.add_argument(
        "--input",
        required=True,
        help="Input mosdepth regions file (gzip-compressed BED format)"
    )
    parser.add_argument(
        "--output",
        required=True,
        help="Output TSV file with normalized CNV values"
    )

    args = parser.parse_args()
    normalize_cnv_depth(args.input, args.output)