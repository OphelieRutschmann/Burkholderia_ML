"""
convert_pangenome_to_matrix.py
------------------------------
Convert Roary pangenome gene presence/absence matrix to feature matrix format.

Usage:
    python convert_pangenome_to_matrix.py \
        --input results/05_Pangenome/gene_presence_absence.csv \
        --output results/genes_presence_absence.tsv
"""

import argparse
import pandas as pd


def convert_roary_matrix(input_csv):
    """
    Convert Roary's gene_presence_absence.csv to feature matrix format.

    Roary output has:
    - First column: gene name
    - Subsequent columns: sample names with gene presence/absence

    Returns:
        - feature_matrix: DataFrame with samples as rows, genes as columns
        - Binary values: 1 = present, 0 = absent
    """

    # Read Roary output
    df = pd.read_csv(input_csv, index_col=0)

    # Transpose: samples as rows, genes as columns
    feature_matrix = df.T

    # Convert to binary (Roary often has gene names for present genes)
    # Replace any non-zero/non-null values with 1, NaN/empty with 0
    feature_matrix = (feature_matrix.notna() & (feature_matrix != '')).astype(int)

    feature_matrix.index.name = "sample"

    print(f"Pangenome matrix: {feature_matrix.shape[0]} samples × {feature_matrix.shape[1]} genes")
    print(f"Core genes (present in all samples): {(feature_matrix.sum(axis=0) == len(feature_matrix)).sum()}")
    print(f"Accessory genes (present in some): {(feature_matrix.sum(axis=0) > 0).sum()}")

    return feature_matrix


def main():
    parser = argparse.ArgumentParser(description="Convert Roary pangenome matrix to feature matrix.")
    parser.add_argument("--input", required=True, help="Path to Roary gene_presence_absence.csv")
    parser.add_argument("--output", required=True, help="Output feature matrix (TSV)")
    args = parser.parse_args()

    print(f"Reading pangenome matrix: {args.input}")
    feature_matrix = convert_roary_matrix(args.input)

    # Save feature matrix
    feature_matrix.to_csv(args.output, sep="\t")
    print(f"\nFeature matrix saved to: {args.output}")
    print(f"Matrix shape: {feature_matrix.shape[0]} samples × {feature_matrix.shape[1]} genes")


if __name__ == "__main__":
    main()