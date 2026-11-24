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

Core SVs are those whose locus is covered (no -1) in at least a given fraction of the
samples (--min-covered-samples, default 0.95). In the few samples where a kept SV's locus
is not covered, the SV is counted as absent (0). Requiring coverage in every sample
removes nearly everything: each strain lacks a different part of the accessory genome.
SVs that are present in all samples or in none carry no information and are dropped.

Usage:
    python convert_sv_to_matrix.py \\
        --vcf results/03_SV/merged_svs.vcf.gz \\
        --samples results/03_SV/merged_svs.samples.txt \\
        --coverage S1=results/03_SV/S1/flanks.thresholds.bed.gz S2=... \\
        --min-frac 0.9 \\
        --min-covered-samples 0.95 \\
        --output results/svs_core.tsv \\
        --metadata results/svs_core_meta.tsv
"""

import argparse
import gzip
import pandas as pd


def read_vcf_records(path):
    """
    Yield (ID, INFO dict) for each record of a (bgzipped) VCF.
    Plain-Python parsing, so no VCF library is needed in the container.
    """
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            info = {}
            for item in fields[7].split(";"):
                key, _, value = item.partition("=")
                info[key] = value if value else True
            yield fields[2], info


def read_flank_coverage(path):
    """
    Read a mosdepth thresholds file (chrom, start, end, region, bases >= depth)
    and return {region_name: fraction of bases at or above the depth threshold}.
    """
    df = pd.read_csv(path, sep="\t")
    length = (df.iloc[:, 2] - df.iloc[:, 1]).clip(lower=1)
    return dict(zip(df.iloc[:, 3], df.iloc[:, 4] / length))


def parse_svlen(value):
    """SV length as a positive integer. SURVIVOR writes long lengths in scientific
    notation (e.g. 1.97908e+06); missing or unreadable values give 0."""
    try:
        return int(abs(float(value)))
    except (TypeError, ValueError):
        return 0


def locus_covered(sv_id, flank_frac, min_frac):
    """Both flanks of the SV (when they exist) must be covered."""
    fracs = [flank_frac[f"{sv_id}__{side}"] for side in ("L", "R") if f"{sv_id}__{side}" in flank_frac]
    return bool(fracs) and all(f >= min_frac for f in fracs)


def extract_core_svs(vcf_path, samples, coverage, min_frac, min_covered_samples):
    """
    Returns:
        - sv_df: DataFrame with samples as rows, core SVs as columns (0/1 presence)
        - sv_meta: dict with SV metadata (svtype, svlen, n_present, frac_samples_covered)
    """
    sv_rows = {}
    sv_meta = {}

    for sv_id, info in read_vcf_records(vcf_path):
        supp_vec = info.get("SUPP_VEC")
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
            "svtype": info.get("SVTYPE", "UNKNOWN"),
            "svlen": parse_svlen(info.get("SVLEN")),
        }

    sv_df_all = pd.DataFrame(sv_rows, index=samples)
    sv_df_all.index.name = "sample"

    # Fraction of samples in which the SV locus is covered (called SVs count as covered)
    frac_covered = (sv_df_all != -1).mean(axis=0)
    n_present = (sv_df_all == 1).sum(axis=0)
    variable = (n_present > 0) & (n_present < len(samples))

    print(f"[SV] Total merged SV events: {sv_df_all.shape[1]}")
    print(f"[SV] Samples: {len(samples)}")

    # Diagnostics: per-sample coverage of the SV loci, and SVs kept at other thresholds
    frac_per_sample = (sv_df_all != -1).mean(axis=1).sort_values()
    print("[SV] Fraction of SV loci covered per sample: "
          f"median {frac_per_sample.median():.3f}, min {frac_per_sample.min():.3f}")
    print("[SV] 10 samples with the lowest coverage of SV loci:")
    for sample, frac in frac_per_sample.head(10).items():
        print(f"[SV]     {sample}\t{frac:.3f}")
    print("[SV] Variable SVs that would be kept, by minimum fraction of samples covered:")
    for t in (1.0, 0.99, 0.95, 0.9, 0.8):
        print(f"[SV]     >= {t:.2f}: {int(((frac_covered >= t) & variable).sum())}")

    covered = frac_covered >= min_covered_samples
    keep = sv_df_all.columns[covered & variable]

    # Uncovered loci in kept SVs are counted as absent
    sv_df = sv_df_all[keep].replace(-1, 0)
    n_filled = int((sv_df_all[keep] == -1).sum().sum())

    sv_meta_core = {
        sv_id: {**sv_meta[sv_id], "n_present": int(n_present[sv_id]),
                "frac_samples_covered": round(float(frac_covered[sv_id]), 4)}
        for sv_id in keep
    }

    print(f"[SV] Removed {(~covered).sum()} SVs whose locus is covered in fewer than "
          f"{min_covered_samples:.0%} of samples")
    print(f"[SV] Removed {(covered & ~variable).sum()} invariant SVs")
    print(f"[SV] Core SV events kept: {len(keep)}")
    print(f"[SV] Uncovered loci counted as absent in kept SVs: {n_filled} "
          f"({n_filled / max(sv_df.size, 1):.2%} of matrix cells)")

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
    parser.add_argument("--min-covered-samples", type=float, default=0.95,
                        help="Minimum fraction of samples in which the SV locus must be covered (default 0.95)")
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
    sv_df, sv_meta = extract_core_svs(args.vcf, samples, coverage, args.min_frac, args.min_covered_samples)

    sv_df.to_csv(args.output, sep="\t")
    print(f"\nCore SV feature matrix saved to: {args.output}")
    print(f"Matrix shape: {sv_df.shape[0]} samples × {sv_df.shape[1]} SVs")

    meta_df = pd.DataFrame(sv_meta, index=["svtype", "svlen", "n_present", "frac_samples_covered"]).T
    meta_df.index.name = "sv_id"
    meta_df.to_csv(args.metadata, sep="\t")
    print(f"SV metadata saved to: {args.metadata}")


if __name__ == "__main__":
    main()
