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
        raise ValueError("No sequences found in {}".format(alignment_path))
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
            snp_id = "snp_core_{}".format(pos + 1)
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
    print("Alignment: {} samples x {} positions".format(len(sequences), seq_length))
    print("Polymorphic sites (SNPs): {}".format(snp_count))
    return feature_df, ref_alleles


def main():
    parser = argparse.ArgumentParser(description="Convert core SNP alignment to feature matrix.")
    parser.add_argument("--alignment", required=True, help="Path to core SNP alignment (FASTA)")
    parser.add_argument("--output", required=True, help="Output feature matrix (TSV)")
    args = parser.parse_args()
    print("Reading alignment: {}".format(args.alignment))
    feature_df, ref_alleles = alignment_to_matrix(args.alignment)
    # Save feature matrix
    feature_df.to_csv(args.output, sep="\t")
    print("\nFeature matrix saved to: {}".format(args.output))
    print("Matrix shape: {} samples x {} SNPs".format(feature_df.shape[0], feature_df.shape[1]))


if __name__ == "__main__":
    main()
