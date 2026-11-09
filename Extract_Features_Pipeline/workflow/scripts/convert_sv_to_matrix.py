"""
convert_sv_to_matrix.py
-----------------------
Build a sample x SV presence/absence matrix from the merged SV VCF.

The merged VCF is produced with `bcftools merge --missing-to-ref`, so an SV
that was not called in a sample is encoded as absent (0). Any remaining
missing genotypes (-1) cause the SV to be dropped. SVs that are present in
all samples or in none carry no information and are dropped as well.

Usage:
    python convert_sv_to_matrix.py \
        --vcf results/03_SV/merged_svs.vcf.gz \
        --output results/svs_core.tsv \
        --metadata results/svs_core_meta.tsv
"""

import argparse
from cyvcf2 import VCF
import pandas as pd


def extract_core_svs(vcf_path):
    """
    Extract informative SVs (complete and variable across samples) from merged VCF.

    Returns:
        - sv_df: DataFrame with samples as rows, SVs as columns (0/1 presence)
        - sv_meta: dict with SV metadata (svtype, svlen)
    """
    vcf = VCF(vcf_path)
    samples = list(vcf.samples)

    sv_rows = {}
    sv_meta = {}

    for variant in vcf:
        svtype = variant.INFO.get("SVTYPE", "UNKNOWN")
        svlen = abs(int(variant.INFO.get("SVLEN", 0) or 0))

        # CHROM_POS_TYPE alone is not unique (e.g. two deletions of different
        # length starting at the same position), so include END and the ALT allele.
        sv_id = f"sv__{variant.CHROM}_{variant.POS}_{variant.end}_{svtype}"
        if sv_id in sv_rows:
            sv_id = f"{sv_id}_{'-'.join(variant.ALT)}"

        presences = []
        for gt in variant.genotypes:
            alleles = gt[:-1]  # last element is the phased flag
            if any(a == -1 for a in alleles):
                presences.append(-1)
            else:
                presences.append(int(any(a > 0 for a in alleles)))

        sv_rows[sv_id] = presences
        sv_meta[sv_id] = {"svtype": svtype, "svlen": svlen}

    vcf.close()

    sv_df_all = pd.DataFrame(sv_rows, index=samples)
    sv_df_all.columns = sv_df_all.columns.astype(str)
    sv_df_all.index.name = "sample"

    # Keep SVs with complete data
    complete = (sv_df_all == -1).sum(axis=0) == 0
    # Keep SVs that vary across samples (not present in all / none)
    n_present = sv_df_all.sum(axis=0)
    variable = (n_present > 0) & (n_present < len(samples))

    keep = sv_df_all.columns[complete & variable]
    sv_df = sv_df_all[keep].copy()
    sv_meta_core = {sv_id: sv_meta[sv_id] for sv_id in keep}

    print(f"[SV] Total SV events: {sv_df_all.shape[1]}")
    print(f"[SV] Removed {(~complete).sum()} SVs with missing data")
    print(f"[SV] Removed {(complete & ~variable).sum()} invariant SVs")
    print(f"[SV] Informative SV events kept: {len(keep)}")
    print(f"[SV] Samples: {len(samples)}")

    return sv_df, sv_meta_core


def main():
    parser = argparse.ArgumentParser(description="Convert merged SV VCF to a presence/absence matrix.")
    parser.add_argument("--vcf", required=True, help="Path to merged SV VCF (merged_svs.vcf.gz)")
    parser.add_argument("--output", required=True, help="Output SV feature matrix (TSV)")
    parser.add_argument("--metadata", required=True, help="Output SV metadata (TSV)")
    args = parser.parse_args()

    print(f"Reading SV VCF: {args.vcf}")
    sv_df, sv_meta = extract_core_svs(args.vcf)

    # Save feature matrix
    sv_df.to_csv(args.output, sep="\t")
    print(f"\nSV feature matrix saved to: {args.output}")
    print(f"Matrix shape: {sv_df.shape[0]} samples × {sv_df.shape[1]} SVs")

    # Save metadata
    meta_df = pd.DataFrame(sv_meta).T
    meta_df.index.name = "sv_id"
    meta_df.to_csv(args.metadata, sep="\t")
    print(f"SV metadata saved to: {args.metadata}")


if __name__ == "__main__":
    main()