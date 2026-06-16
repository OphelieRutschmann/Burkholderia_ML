# Copy number variation Feature Extraction — `variant_calling_cnv_depth.smk`

## Overview

This Snakemake module computes **copy-number variation (CNV) features** per genomic region for every sample that has passed upstream QC. 
The final output is a single tab-separated feature matrix with one row per sample and one column per genomic region.

---

## Rules

### 1. `prokka` (Reference genome annotation)

| Item | Value |
|---|---|
| Input | `ref_genome.fasta` |
| Output | `prokka.gff`, `prokka.gbk`, `prokka.faa`, `prokka.ffn` |
| Container | `prokka.sif` |

Annotates the reference genome using **Prokka**. Produces a GFF file and other annotation files. The GFF produced here is used by the downstream rule "prepare_regions" to build a target BED file of genomic regions to be analyzed.
---

### 2. `prepare_regions`: Build the target BED file

| Item | Value |
|---|---|
| Input | `prokka/prokka.gff` (from `prokka` rule), `ref_genome.fasta.fai` |
| Output | `regions.bed` |
| Container | `bedtools.sif` |
| Config key | `cnv_depth.promoter_window` (default: 200 bp) |

Constructs a BED file of all genomic regions to be depth-profiled. Each region label has the format:
```
<ID>__<Name>__<feature_type>
```

Two passes are made over the GFF:

**Pass 1: Extract all features.** Every annotated feature (gene, CDS, rRNA, tRNA, etc.) is extracted as a BED interval. Labels encode the GFF `ID`, `Name`, and feature type (gene, CDS, rRNA, tRNA, etc.).

**Pass 2: Extract promoter windows.** For CDS and gene features only, `bedtools flank` creates a new interval upstream of the region of interest by `promoter_window` bp (default 200bp). This corresponds approximately to the promoter region. These new intervals receive the suffix `__promoter` in their label.

Both sets are concatenated, sorted with `bedtools sort`, and deduplicated with `uniq`.
---

### 3. `mosdepth_regions`: Per-sample region depth calculation

| Item | Value |
|---|---|
| Input | `{sample}.bam`, `{sample}.bam.bai`, `regions.bed` |
| Output | `mosdepth.regions.bed.gz` |
| Container | `mosdepth.sif` |

Runs **mosdepth** in region mode (`--by`) over the BED file for each sample, computing mean depth per region. Per-base depth output is suppressed (`--no-per-base`) for efficiency.
---

### 4. `normalize_region_depth`: Compute log₂ CNV ratio

| Item | Value |
|---|---|
| Input | `mosdepth.regions.bed.gz` |
| Output | `cnv_per_region.tsv` |
| Runtime | Python (`pandas`, `numpy`) |

Normalises raw depth values to produce a CNV estimate per region:

1. **Baseline estimation.** The per-chromosome median depth is computed using only `__CDS` and `__gene` features, which represent the core coding content and provide a stable, unbiased baseline. Promoter regions are intentionally excluded because they may carry systematically different coverage.

2. **CNV ratio.** Each region's depth is divided by its chromosome median:

   ```
   cnv_ratio = depth / chrom_median
   ```

3. **Log₂ transformation.** The ratio is log₂-transformed for symmetry around 0 (diploid baseline = 0, duplication > 0, deletion < 0). Zero-depth regions are replaced with `NaN` before the log to avoid `-inf` values in the feature matrix.

Output columns:

| Column | Description |
|---|---|
| `region` | Region label (`ID__Name__type`) |
| `chrom` | Chromosome / contig |
| `depth` | Raw mean depth from mosdepth |
| `chrom_median` | Chromosome baseline (CDS+gene median) |
| `cnv_ratio` | `depth / chrom_median` |
| `log2cnv` | `log2(cnv_ratio)` — the ML feature value |
---

### 5. `merge_depth_cnv`: Assemble the feature matrix

| Item | Value |
|---|---|
| Input | All `cnv_per_region.tsv` files (QC-passed samples only) |
| Output | `cnv_depth_matrix.tsv` |
| Runtime | Python (`pandas`) |

Merges per-sample `log2cnv` vectors into a single matrix:

- Each sample contributes one column of `log2cnv` values indexed by region.
- Columns are concatenated, then transposed so the final shape is **(samples × regions)**.
- All column names are prefixed with `cnv_depth__` to namespace the features when the matrix is joined with other feature sets.
- The row index is named `sample`.
---

## Output Files

```
results/04_CNV_depth/
├── {sample}/
│   ├── mosdepth.regions.bed.gz       # raw per-region depth (mosdepth output)
│   └── cnv_per_region.tsv            # normalised CNV ratios per region
└── cnv_depth_matrix.tsv              # final ML feature matrix (samples × regions)
```

`cnv_depth_matrix.tsv` schema:

| Column | Type | Description |
|---|---|---|
| `sample` (index) | string | Sample identifier |
| `cnv_depth__<ID>__<Name>__<type>` | float | log₂ CNV ratio for that region (`NaN` if depth = 0) |

---

## Configuration

All keys can be found in `config`:

```yaml
resources:
  general:
    cpus: 8
    mem_mb: 16000
    runtime: 120

cnv_depth:
  promoter_window: 200  # bp upstream of CDS/gene to treat as promoter
```

`get_passed_samples()` is imported from `../rules/common.smk` and returns the list of samples that have passed upstream QC filters.
---

## Dependencies

| Tool | Purpose | Provided via |
|---|---|---|
| Prokka | Reference genome annotation | `prokka.sif` |
| bedtools | Interval manipulation, flank | `bedtools.sif` |
| mosdepth | Region-level depth profiling | `mosdepth.sif` |
| pandas | Data wrangling in Python rules | environment |
| numpy | `log2`, `NaN` handling | environment |