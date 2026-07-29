"""
convert_snp_alignment_to_matrix.py
------------------------------
Convert a core SNP alignment (FASTA) to a feature matrix (TSV).

Usage:
    python convert_snp_alignment_to_matrix.py \
        --alignment results/04_Phylogeny/core_SNP/clean.full.aln \
        --output results/snps_core.tsv
"""

import argparse
from Bio import SeqIO
import pandas as pd
import numpy as np


def alignment_to_matrix(alignment_path):
    """
    Convert FASTA alignment to SNP feature matrix.

    Returns:
        - feature_matrix: DataFrame with samples as rows, SNP positions as columns
        - ref_alleles: dict mapping position to reference allele
    """
    sequences = {}
    ref_allele = None

    # Read FASTA alignment
    for record in SeqIO.parse(alignment_path, "fasta"):
        sample_id = record.id
        seq = str(record.seq).upper()
        sequences[sample_id] = seq

    if not sequences:
        raise ValueError(f"No sequences found in {alignment_path}")

    # Get sequence length
    seq_length = len(list(sequences.values())[0])

    # Extract SNPs (positions with variation)
    matrix = {}
    ref_alleles = {}
    snp_count = 0

    for pos in range(seq_length):
        alleles = set()
        for sample_id, seq in sequences.items():
            allele = seq[pos]
            if allele != '-':  # Skip gaps
                alleles.add(allele)

        # Keep positions with variation (polymorphic sites)
        if len(alleles) > 1:
            snp_count += 1
            snp_id = f"snp_core_{pos+1}"
            ref_alleles[snp_id] = list(alleles)[0]  # First allele as reference

            # Create binary genotypes (0 = ref, 1 = alt)
            matrix[snp_id] = []
            for sample_id in sequences.keys():
                allele = sequences[sample_id][pos]
                if allele == '-':
                    matrix[snp_id].append(-1)  # Missing
                elif allele == ref_alleles[snp_id]:
                    matrix[snp_id].append(0)
                else:
                    matrix[snp_id].append(1)

    # Create DataFrame
    feature_df = pd.DataFrame(matrix, index=list(sequences.keys()))
    feature_df.index.name = "sample"

    print(f"Alignment: {len(sequences)} samples × {seq_length} positions")
    print(f"Polymorphic sites (SNPs): {snp_count}")

    return feature_df, ref_alleles


def main():
    parser = argparse.ArgumentParser(description="Convert core SNP alignment to feature matrix.")
    parser.add_argument("--alignment", required=True, help="Path to core SNP alignment (FASTA)")
    parser.add_argument("--output", required=True, help="Output feature matrix (TSV)")
    args = parser.parse_args()

    print(f"Reading alignment: {args.alignment}")
    feature_df, ref_alleles = alignment_to_matrix(args.alignment)

    # Save feature matrix
    feature_df.to_csv(args.output, sep="\t")
    print(f"\nFeature matrix saved to: {args.output}")
    print(f"Matrix shape: {feature_df.shape[0]} samples × {feature_df.shape[1]} SNPs")


if __name__ == "__main__":
    main()