"""
create_feature_matrix.py
------------------------
Combines SNP, SV, and CNV features into a single per-sample feature matrix for machine learning.

Features extracted:
- SNPs: Individual variants with HIGH/MODERATE/MODIFIER impact (excludes LOW)
- SVs: Structural variants (TD, INV, LI, D types) - presence/absence and counts
- CNVs: Depth-based copy-number variations per genomic region

Usage:
    python create_feature_matrix.py \
        --snp  results/02_Variant_SNPs/snps.ann.vcf \
        --sv   results/03_SV/merged_svs.vcf.gz \
        --cnv-depth  results/04_CNV/cnv_depth_matrix.tsv \
        --out  results/feature_matrix.tsv

Outputs:
    - <out> feature matrix (samples vs features)
    - <out>_snp_meta.tsv SNP variant metadata (var_id vs gene/impact/effect)
    - <out>_sv_meta.tsv SV metadata (sv_id vs svtype/length)
"""

import argparse
from collections import defaultdict

import pandas as pd
from cyvcf2 import VCF

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

ANN_FIELDS = [
    "allele", "effect", "putative_impact", "gene_name", "gene_id",
    "feature_type", "feature_id", "transcript_biotype", "rank",
    "hgvs_c", "hgvs_p", "cdna_pos", "cds_pos", "protein_pos", "distance", "info"
]

def parse_ann(info_ann: str) -> list[dict]:
    """Parse the ANN field from the SNPEff-annotated VCF into a list of dicts."""
    annotations = []
    for entry in info_ann.split(","):
        parts = entry.split("|")
        ann = dict(zip(ANN_FIELDS, parts + [""] * (len(ANN_FIELDS) - len(parts))))
        annotations.append(ann)
    return annotations

# ---------------------------------------------------------------------------
# SNP feature extraction
# ---------------------------------------------------------------------------

def extract_snp_genotypes(vcf_path: str):
    """
    Description:
    Creates a genotype presence/absence matrix and a feature metadata dictionary
    from a multi-sample VCF annotated with SNPEff.

    Outputs:
    1) geno_df: DataFrame with samples as rows and variants (var_id) as columns.
       Values:
           -1  genotype missing for this sample at this site
            0  sample has the reference allele (or a different alt)
            1  sample has this specific alt allele

    2) feature_meta: dict keyed by var_id containing gene, impact, and effect
       for each variant (strongest SNPEff annotation across transcripts).
    """
    vcf = VCF(vcf_path)
    samples = list(vcf.samples)

    geno_rows = {}      # var_id -> list of presences, one per sample
    feature_meta = {}   # var_id -> {gene, impact, effect}

    for variant in vcf:
        chrom = variant.CHROM
        pos = variant.POS
        ref = variant.REF
        alts = variant.ALT
        genotypes = variant.genotypes

        a1_calls = [gt[0] for gt in genotypes]  # haploid: first allele only

        # Parse SNPEff annotations and index them by alt allele
        ann_str = variant.INFO.get("ANN")
        annotations = parse_ann(ann_str) if ann_str else []

        ann_by_alt = {}
        for ann in annotations:
            alt_allele = ann["allele"]
            ann_by_alt.setdefault(alt_allele, []).append(ann)

        # One feature per alt allele
        for alt_idx, alt in enumerate(alts):
            target = alt_idx + 1
            var_id = f"snp__{chrom}_{pos}_{ref}_{alt}"

            # 1. Genotype matrix
            presences = [
                -1 if a1 == -1 else int(a1 == target)
                for a1 in a1_calls
            ]
            geno_rows[var_id] = presences

            # 2. Metadata: pick strongest SNPEff annotation for this allele
            feature_meta[var_id] = {
                "gene": None,
                "impact": None,
                "effect": None
            }

            ann_list = ann_by_alt.get(alt, [])
            if ann_list:
                priority = {"HIGH": 4, "MODERATE": 3, "LOW": 2, "MODIFIER": 1}
                best = max(
                    ann_list,
                    key=lambda x: priority.get(x.get("putative_impact", ""), 0)
                )
                feature_meta[var_id] = {
                    "gene":   best.get("gene_name", ""),
                    "impact": best.get("putative_impact", ""),
                    "effect": best.get("effect", "")
                }

    vcf.close()

    geno_df = pd.DataFrame(geno_rows, index=samples)
    geno_df.index.name = "sample"

    print(f"[SNP] {len(geno_rows)} variant sites across {len(samples)} samples")
    return geno_df, feature_meta

# ---------------------------------------------------------------------------
# CNV feature extraction
# ---------------------------------------------------------------------------

SVTYPES = ["TD", "INV", "LI", "D"]

def extract_cnv_features(cnv_vcf_path: str) -> pd.DataFrame:
    """
    Returns a DataFrame (samples X CNV features) with:
      - Binary presence/absence per SV event   (cnv__<CHROM>_<POS>_<SVTYPE>)
      - Per-sample, per-SVTYPE count           (cnv_count_<SVTYPE>)
      - Per-sample, per-SVTYPE total bp        (cnv_bp_<SVTYPE>)
    """
    vcf = VCF(cnv_vcf_path)
    samples = vcf.samples

    sv_rows = {}
    sv_counts = {s: defaultdict(int) for s in samples}
    sv_bp = {s: defaultdict(int) for s in samples}

    for variant in vcf:
        chrom, pos = variant.CHROM, variant.POS
        svtype = variant.INFO.get("SVTYPE", "UNKNOWN")
        svlen  = abs(int(variant.INFO.get("SVLEN", 0)))

        sv_id = f"cnv__{chrom}_{pos}_{svtype}"
        presences = []

        for gt in variant.genotypes:
            a1, a2 = gt[0], gt[1]
            present = int(a1 > 0 or a2 > 0) if (a1 != -1 and a2 != -1) else -1
            presences.append(present)

        sv_rows[sv_id] = presences

        for sample_idx, sample in enumerate(samples):
            if presences[sample_idx] > 0:
                sv_counts[sample][f"cnv_count_{svtype}"] += 1
                sv_bp[sample][f"cnv_bp_{svtype}"]        += svlen

    vcf.close()

    sv_df = pd.DataFrame(sv_rows, index=samples)
    sv_df.columns = sv_df.columns.astype(str)

    count_df = pd.DataFrame(sv_counts).T.fillna(0).astype(int)
    bp_df    = pd.DataFrame(sv_bp).T.fillna(0).astype(int)

    # Ensure all SV types are present even if none were observed
    for svtype in SVTYPES:
        for df, prefix in [(count_df, "cnv_count_"), (bp_df, "cnv_bp_")]:
            col = f"{prefix}{svtype}"
            if col not in df.columns:
                df[col] = 0

    cnv_features = pd.concat([sv_df, count_df, bp_df], axis=1)
    cnv_features.index.name = "sample"
    print(f"[CNV] {len(sv_rows)} SV events | "
          f"{count_df.shape[1]} count features | "
          f"{bp_df.shape[1]} bp-burden features")
    return cnv_features


# ---------------------------------------------------------------------------
# Structural Variant feature extraction
# ---------------------------------------------------------------------------

def extract_sv_features(sv_vcf_path: str) -> tuple[pd.DataFrame, dict]:
    """
    Extract structural variant features from a merged SV VCF.

    Returns:
    1) sv_df: DataFrame with samples as rows and SV variants as columns.
       Values: -1 (missing), 0 (absent), 1 (present)
    2) sv_meta: dict keyed by sv_id containing svtype and svlen
    """
    vcf = VCF(sv_vcf_path)
    samples = list(vcf.samples)

    sv_rows = {}
    sv_meta = {}
    sv_types = ["TD", "INV", "LI", "D"]

    for variant in vcf:
        chrom = variant.CHROM
        pos = variant.POS
        svtype = variant.INFO.get("SVTYPE", "UNKNOWN")
        svlen = abs(int(variant.INFO.get("SVLEN", 0)))

        sv_id = f"sv__{chrom}_{pos}_{svtype}"
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

    sv_df = pd.DataFrame(sv_rows, index=samples)
    sv_df.columns = sv_df.columns.astype(str)
    sv_df.index.name = "sample"

    print(f"[SV] {len(sv_rows)} SV events across {len(samples)} samples")
    return sv_df, sv_meta


# ---------------------------------------------------------------------------
# Depth-based CNV matrix loading
# ---------------------------------------------------------------------------

def load_cnv_depth_matrix(cnv_depth_path: str) -> pd.DataFrame:
    """Load the pre-computed depth-based CNV feature matrix."""
    cnv_df = pd.read_csv(cnv_depth_path, sep="\t", index_col=0)
    cnv_df.index.name = "sample"
    print(f"[CNV-depth] Loaded {cnv_df.shape[1]} depth-based CNV features")
    return cnv_df


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description="Build ML feature matrix from SNP, SV, and CNV features.")
    parser.add_argument("--snp", required=True, help="Path to annotated SNP VCF (snps.ann.vcf)")
    parser.add_argument("--sv", required=True, help="Path to merged SV VCF (merged_svs.vcf.gz)")
    parser.add_argument("--cnv-depth", required=True, help="Path to depth-based CNV matrix (cnv_depth_matrix.tsv)")
    parser.add_argument("--out", required=True, help="Output TSV path for the feature matrix")
    args = parser.parse_args()

    print("Extracting SNP features...")
    snp_df, snp_meta = extract_snp_genotypes(args.snp)

    # Drop variants whose strongest effect is LOW
    keep = [
        var_id for var_id, meta in snp_meta.items()
        if meta["impact"] in ("HIGH", "MODERATE", "MODIFIER")
    ]
    snp_df = snp_df[keep]
    snp_meta = {var_id: snp_meta[var_id] for var_id in keep}
    print(f"[SNP] {len(keep)} variants kept after dropping LOW")

    print("Extracting SV features...")
    sv_df, sv_meta = extract_sv_features(args.sv)

    print("Loading depth-based CNV features...")
    cnv_depth_df = load_cnv_depth_matrix(args.cnv_depth)

    # Align on shared samples across all three feature types
    shared = snp_df.index.intersection(sv_df.index).intersection(cnv_depth_df.index)
    all_samples = snp_df.index.union(sv_df.index).union(cnv_depth_df.index)

    missing_snp = all_samples.difference(snp_df.index).tolist()
    missing_sv = all_samples.difference(sv_df.index).tolist()
    missing_cnv = all_samples.difference(cnv_depth_df.index).tolist()

    if missing_snp:
        print(f"[WARN] {len(missing_snp)} samples missing from SNP VCF")
    if missing_sv:
        print(f"[WARN] {len(missing_sv)} samples missing from SV VCF")
    if missing_cnv:
        print(f"[WARN] {len(missing_cnv)} samples missing from CNV depth matrix")

    if len(shared) == 0:
        raise ValueError("No shared samples found across all feature files!")

    print(f"[INFO] Using {len(shared)} shared samples")

    # Combine features
    matrix = pd.concat([
        snp_df.loc[shared],
        sv_df.loc[shared],
        cnv_depth_df.loc[shared]
    ], axis=1)
    matrix.index.name = "sample"
    matrix.to_csv(args.out, sep="\t")

    # Save metadata alongside the matrix
    snp_meta_path = args.out.replace(".tsv", "_snp_meta.tsv")
    pd.DataFrame(snp_meta).T.to_csv(snp_meta_path, sep="\t", index_label="var_id")

    sv_meta_path = args.out.replace(".tsv", "_sv_meta.tsv")
    pd.DataFrame(sv_meta).T.to_csv(sv_meta_path, sep="\t", index_label="sv_id")

    print(f"\nFeature matrix written to: {args.out}")
    print(f"SNP metadata written to: {snp_meta_path}")
    print(f"SV metadata written to: {sv_meta_path}")
    print(f"\nMatrix shape: {matrix.shape[0]} samples × {matrix.shape[1]} features")
    print(f"Feature breakdown:")
    print(f"SNP features: {snp_df.shape[1]}")
    print(f"SV features: {sv_df.shape[1]}")
    print(f"CNV (depth) features: {cnv_depth_df.shape[1]}")


if __name__ == "__main__":
    main()