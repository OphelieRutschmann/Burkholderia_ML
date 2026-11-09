"""
convert_sv_to_matrix.py
-----------------------
Build a sample x SV presence/absence matrix of core SVs from the SURVIVOR-merged VCF.

For every merged SV and every sample:
  1  the SV was called in the sample (INFO/SUPP_VEC)
  0  the SV was not called, and both flanks of the SV locus are covered by reads
     in that sample (the locus exists in the strain, so the SV is truly absent)
 -1  not called and a flank is not covered: the locus may be missing from that
     strain's genome (accessory region), so absence cannot be told apart from
     "region not present"

Core SVs are those without -1 in any sample. SVs that are present in all samples
or in none carry no information and are dropped as well.

Usage:
    python convert_sv_to_matrix.py \\
        --vcf results/03_SV/merged_svs.vcf.gz \\
        --samples results/03_SV/merged_svs.samples.txt \\
        --coverage S1=results/03_SV/S1/flanks.thresholds.bed.gz S2=... \\
        --min-frac 0.9 \\
        --output results/svs_core.tsv \\
        --metadata results/svs_core_meta.tsv
"""

import argparse
from cyvcf2 import VCF
import pandas as pd


def read_flank_coverage(path):
    """
    Read a mosdepth thresholds file (chrom, start, end, region, bases >= depth)
    and return {region_name: fraction of bases at or above the depth threshold}.
    """
    df = pd.read_csv(path, sep="\t")
    length = (df.iloc[:, 2] - df.iloc[:, 1]).clip(lower=1)
    return dict(zip(df.iloc[:, 3], df.iloc[:, 4] / length))


def locus_covered(sv_id, flank_frac, min_frac):
    """Both flanks of the SV (when they exist) must be covered."""
    fracs = [flank_frac[f"{sv_id}__{side}"] for side in ("L", "R") if f"{sv_id}__{side}" in flank_frac]
    return bool(fracs) and all(f >= min_frac for f in fracs)


def extract_core_svs(vcf_path, samples, coverage, min_frac):
    """
    Returns:
        - sv_df: DataFrame with samples as rows, core SVs as columns (0/1 presence)
        - sv_meta: dict with SV metadata (svtype, svlen, n_present)
    """
    sv_rows = {}
    sv_meta = {}

    vcf = VCF(vcf_path)
    for variant in vcf:
        sv_id = variant.ID
        supp_vec = variant.INFO.get("SUPP_VEC")
        if supp_vec is None or len(supp_vec) != len(samples):
            raise ValueError(
                f"{sv_id}: SUPP_VEC '{supp_vec}' does not match the {len(samples)} samples"
            )

        values = []
        for sample, called in zip(samples, supp_vec):
            if called == "1":
                values.append(1)
            elif locus_covered(sv_id, coverage[sample], min_frac):
                values.append(0)
            else:
                values.append(-1)

        sv_rows[sv_id] = values
        sv_meta[sv_id] = {
            "svtype": variant.INFO.get("SVTYPE", "UNKNOWN"),
            "svlen": abs(int(variant.INFO.get("SVLEN", 0) or 0)),
        }
    vcf.close()

    sv_df_all = pd.DataFrame(sv_rows, index=samples)
    sv_df_all.index.name = "sample"

    complete = (sv_df_all == -1).sum(axis=0) == 0
    n_present = (sv_df_all == 1).sum(axis=0)
    variable = (n_present > 0) & (n_present < len(samples))

    keep = sv_df_all.columns[complete & variable]
    sv_df = sv_df_all[keep].copy()
    sv_meta_core = {sv_id: {**sv_meta[sv_id], "n_present": int(n_present[sv_id])} for sv_id in keep}

    print(f"[SV] Total merged SV events: {sv_df_all.shape[1]}")
    print(f"[SV] Removed {(~complete).sum()} SVs whose locus is not covered in every sample")
    print(f"[SV] Removed {(complete & ~variable).sum()} invariant SVs")
    print(f"[SV] Core SV events kept: {len(keep)}")
    print(f"[SV] Samples: {len(samples)}")

    return sv_df, sv_meta_core


def main():
    parser = argparse.ArgumentParser(description="Convert merged SV VCF to a core SV presence/absence matrix.")
    parser.add_argument("--vcf", required=True, help="SURVIVOR-merged SV VCF (merged_svs.vcf.gz)")
    parser.add_argument("--samples", required=True,
                        help="Sample names, one per line, in the order of the SURVIVOR input list (SUPP_VEC order)")
    parser.add_argument("--coverage", required=True, nargs="+",
                        help="SAMPLE=path pairs to mosdepth flank thresholds files")
    parser.add_argument("--min-frac", type=float, default=0.9,
                        help="Fraction of flank bases that must reach the depth threshold (default 0.9)")
    parser.add_argument("--output", required=True, help="Output SV feature matrix (TSV)")
    parser.add_argument("--metadata", required=True, help="Output SV metadata (TSV)")
    args = parser.parse_args()

    with open(args.samples) as f:
        samples = [line.strip() for line in f if line.strip()]

    coverage = {}
    for item in args.coverage:
        sample, path = item.split("=", 1)
        coverage[sample] = read_flank_coverage(path)
    missing = set(samples) - set(coverage)
    if missing:
        raise ValueError(f"No coverage file given for samples: {sorted(missing)}")

    print(f"Reading SV VCF: {args.vcf}")
    sv_df, sv_meta = extract_core_svs(args.vcf, samples, coverage, args.min_frac)

    sv_df.to_csv(args.output, sep="\t")
    print(f"\nCore SV feature matrix saved to: {args.output}")
    print(f"Matrix shape: {sv_df.shape[0]} samples × {sv_df.shape[1]} SVs")

    meta_df = pd.DataFrame(sv_meta, index=["svtype", "svlen", "n_present"]).T
    meta_df.index.name = "sv_id"
    meta_df.to_csv(args.metadata, sep="\t")
    print(f"SV metadata saved to: {args.metadata}")


if __name__ == "__main__":
    main()
