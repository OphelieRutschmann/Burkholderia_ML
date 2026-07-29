"""
extract_core_svs.py
-------------------
Extract core SVs (present in most strains) from merged SV VCF.

Core SVs are those with complete data (non-missing) across all samples,
avoiding confounding from variable presence of genomic regions.

Usage:
    python extract_core_svs.py \
        --vcf results/03_SV/merged_svs.vcf.gz \
        --output results/svs_core.tsv \
        --metadata results/svs_core_meta.tsv
"""

import argparse
from cyvcf2 import VCF
import pandas as pd


def extract_core_svs(vcf_path):
    """
    Extract core SVs (with complete data across all samples) from merged VCF.

    Returns:
        - sv_df: DataFrame with samples as rows, SVs as columns (0/1 genotypes)
        - sv_meta: dict with SV metadata (svtype, svlen)
    """
    vcf = VCF(vcf_path)
    samples = list(vcf.samples)

    sv_rows = {}
    sv_meta = {}
    all_svs = []

    for variant in vcf:
        chrom = variant.CHROM
        pos = variant.POS
        svtype = variant.INFO.get("SVTYPE", "UNKNOWN")
        svlen = abs(int(variant.INFO.get("SVLEN", 0)))

        sv_id = f"sv__{chrom}_{pos}_{svtype}"
        all_svs.append(sv_id)

        presences = []
        for gt in variant.genotypes:
            a1, a2 = gt[0], gt[1]
            present = int(a1 > 0 or a2 > 0) if (a1 != -1 and a2 != -1) else -1
            presences.append(present)

        sv_rows[sv_id] = presences
        sv_meta[sv_id] = {
            "svtype": svtype,
            "svlen": svlen
        }

    vcf.close()

    # Create DataFrame with all SVs
    sv_df_all = pd.DataFrame(sv_rows, index=samples)
    sv_df_all.columns = sv_df_all.columns.astype(str)
    sv_df_all.index.name = "sample"

    # Filter to core SVs (complete data, no missing values)
    missing_per_sv = (sv_df_all == -1).sum(axis=0)
    core_svs = sv_df_all.columns[missing_per_sv == 0]

    sv_df = sv_df_all[core_svs].copy()
    sv_meta_core = {sv_id: sv_meta[sv_id] for sv_id in core_svs if sv_id in sv_meta}

    print(f"[SV] Total SV events: {len(all_svs)}")
    print(f"[SV] Core SV events (complete data): {len(core_svs)}")
    print(f"[SV] Removed {len(all_svs) - len(core_svs)} SVs with missing data")
    print(f"[SV] Samples: {len(samples)}")

    return sv_df, sv_meta_core


def main():
    parser = argparse.ArgumentParser(description="Extract core SVs from merged VCF.")
    parser.add_argument("--vcf", required=True, help="Path to merged SV VCF (merged_svs.vcf.gz)")
    parser.add_argument("--output", required=True, help="Output SV feature matrix (TSV)")
    parser.add_argument("--metadata", required=True, help="Output SV metadata (TSV)")
    args = parser.parse_args()

    print(f"Reading SV VCF: {args.vcf}")
    sv_df, sv_meta = extract_core_svs(args.vcf)

    # Save feature matrix
    sv_df.to_csv(args.output, sep="\t")
    print(f"\nCore SV feature matrix saved to: {args.output}")
    print(f"Matrix shape: {sv_df.shape[0]} samples × {sv_df.shape[1]} SVs")

    # Save metadata
    meta_df = pd.DataFrame(sv_meta).T
    meta_df.index.name = "sv_id"
    meta_df.to_csv(args.metadata, sep="\t")
    print(f"SV metadata saved to: {args.metadata}")


if __name__ == "__main__":
    main()