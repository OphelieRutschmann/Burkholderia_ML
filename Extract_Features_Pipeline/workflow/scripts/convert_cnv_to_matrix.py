"""
convert_cnv_to_matrix.py
--------------------------------
Convert CNV depth ratios to gene presence/absence with copy number calls.

Input: Individual CNV per-region TSV files from mosdepth
Output: Binary/ternary gene presence/absence matrix

Copy number calling:
  0 = gene absent (log2cnv < -1.0)
  1 = gene present, normal copy (log2cnv in [-1.0, 0.5])
  2+ = multiple copies (log2cnv > 0.5)

Usage:
    python convert_cnv_depth_to_presence.py \
        --samples sample1 sample2 sample3 \
        --input-dir results/04_CNV \
        --output results/cnvs_gene_presence.tsv
"""

import argparse
import pandas as pd
import numpy as np


def convert_cnv_files(input_dir, samples):
    """
    Convert individual CNV TSV files to copy number matrix.

    Reads cnv_per_region.tsv files from each sample and converts
    log2cnv values to copy number calls.
    """
    copy_dfs = []

    for sample in samples:
        tsv_path = f"{input_dir}/{sample}/cnv_per_region.tsv"

        try:
            df = pd.read_csv(tsv_path, sep="\t")
        except FileNotFoundError:
            print(f"WARNING: File not found: {tsv_path}")
            continue

        # Convert log2cnv to copy number
        # 0 = gene absent (log2 < -1.0)
        # 1 = gene present, normal copy (log2 in [-1.0, 0.5])
        # 2+ = multiple copies (log2 > 0.5)
        df["copy_number"] = pd.cut(
            df["log2cnv"],
            bins=[-np.inf, -1.0, 0.5, np.inf],
            labels=[0, 1, 2],
            right=False
        ).astype(int)

        df_copy = df[["region", "copy_number"]].copy()
        df_copy.columns = ["region", sample]
        df_copy = df_copy.set_index("region")
        copy_dfs.append(df_copy)

    # Concatenate all samples
    copy_matrix = pd.concat(copy_dfs, axis=1).T
    copy_matrix.index.name = "sample"
    copy_matrix.columns = [f"gene__{c}" for c in copy_matrix.columns]

    return copy_matrix


def main():
    parser = argparse.ArgumentParser(description="Convert CNV depth ratios to gene presence/absence.")
    parser.add_argument("--samples", required=True, nargs="+", help="List of sample IDs")
    parser.add_argument("--input-dir", required=True, help="Path to CNV results directory")
    parser.add_argument("--output", required=True, help="Output gene presence/absence matrix (TSV)")
    args = parser.parse_args()

    print(f"Converting CNV depth ratios to gene presence/absence...")
    copy_matrix = convert_cnv_files(args.input_dir, args.samples)

    # Save
    copy_matrix.to_csv(args.output, sep="\t")
    print(f"\nGene presence/absence matrix saved to: {args.output}")
    print(f"Matrix shape: {copy_matrix.shape[0]} samples × {copy_matrix.shape[1]} genes")

    # Summary
    print(f"\nCopy number distribution:")
    for col in copy_matrix.columns[:5]:  # Show first 5 genes
        dist = copy_matrix[col].value_counts().sort_index()
        print(f"  {col}: {dict(dist)}")


if __name__ == "__main__":
    main()