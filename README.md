# Burkholderia_ML

Predicting co-trimoxazole (trimethoprim–sulfamethoxazole, **SXT**) resistance in *Burkholderia pseudomallei* from whole-genome sequencing data using machine learning.

The repository has two parts:

1. **`Extract_Features_Pipeline/`** – a Snakemake workflow that turns paired-end Illumina reads into genomic feature matrices (core-genome SNPs, gene/promoter copy-number depth, optionally structural variants) plus a core-SNP phylogeny.
2. **`snp_ML/` and `cnv_ML/`** – Jupyter notebooks that clean those matrices, select features and train/compare classifiers (logistic regression, KNN, SVM, decision tree, XGBoost), with final performance estimated by sequence-type–aware nested cross-validation.

---

## Dataset

- **227 clinical *B. pseudomallei* isolates** (all from patients, mostly northeastern Thailand – Khon Kaen, Nong Khai, Sisaket, …), identified by SRA accession.
- Phenotype: binary SXT label (`SXT_binary`) parsed from the `resistance_profile` field – **123 resistant / 104 susceptible**.
- 59 distinct sequence types (STs). STs are used to build train/test splits so that test isolates come from lineages not seen during training, limiting phylogenetic leakage.

| File | Content |
|---|---|
| `genome_metadata.csv` | Raw metadata (accession, source, species, resistance profile, ST, location) |
| `genome_metadata_cleaned.csv` | Filtered metadata for the 227 isolates used |
| `itol_SXT` | iTOL binary dataset to annotate the phylogeny with SXT resistance |

---

## Repository structure

```
Burkholderia_ML/
├── Extract_Features_Pipeline/      # Snakemake feature-extraction workflow
│   ├── Snakefile
│   ├── run_pipeline.sh             # launches the workflow on SLURM
│   ├── create_samplesheet.sh       # builds samples.tsv from a reads folder
│   ├── pull_containers.sh          # pulls/builds Singularity images
│   ├── Burkholderia_ML.yml         # conda env (Snakemake 7.32 + Apptainer)
│   └── workflow/
│       ├── modules/                # QC, databases, SNPs, SVs, CNVs, phylogeny
│       ├── rules/common.smk        # shared helper functions
│       ├── scripts/                # Python helpers (matrix conversion, normalisation…)
│       └── containers/python_env.def
├── snp_ML/                         # ML on core-genome SNP features
│   ├── Data/                       # snps.tsv, snps_cleaned.csv, SXT_cleaned.csv …
│   ├── Figures/
│   └── 01_… – 09_….ipynb
├── cnv_ML/                         # ML on copy-number (read-depth) features
│   ├── Data/                       # cnv_depth_matrix.tsv, train/test splits, clusters
│   ├── Figures/                    # ROC, confusion matrices, biomarkers, annotated trees
│   ├── Results/                    # model metrics, predictions, biomarker rankings
│   ├── nested_cv_results.csv
│   └── 01_… – 10_….ipynb
├── Benchmarking/AbR_reports/       # per-isolate rule-based AMR prediction reports
├── genome_metadata*.csv
└── itol_SXT
```

---

## 1. Feature extraction pipeline

### Workflow overview

| Step | Module | Tools |
|---|---|---|
| Read QC & trimming | `QC_illumina.smk` | FastQC, fastp, MultiQC |
| Mapping & coverage | `QC_illumina.smk` | BWA, SAMtools |
| Sample filtering (checkpoint) | `filter_reads.py` | Q30 rate > 0.90, GC content 0.67 ± 0.03, > 75 % reads surviving trimming, minimum coverage |
| Reference annotation | `databases.smk`, `variant_calling_snps.smk` | Bakta (light DB) |
| SNP calling | `variant_calling_snps.smk` | Snippy (`--mincov 10 --minfrac 0.9`), bcftools |
| Core-SNP alignment & matrix | `variant_calling_snps.smk` | snippy-core, snpEff annotation |
| CNV (read depth) | `variant_calling_cnvs.smk` | bedtools, mosdepth; depth per gene/CDS and its upstream promoter window, normalised to log2 ratio vs. chromosome median |
| Structural variants *(optional, disabled by default)* | `variant_calling_svs.smk` | Picard, Pindel |
| Phylogeny | `phylogeny.smk` | IQ-TREE 2 (GTR+G) on the core alignment |

All tools run in Singularity/Apptainer containers stored in `workflow/containers/`.

### Main outputs

```
results/02_Variant_SNPs/core_SNP/core_snps.ann.vcf    # annotated core SNPs
results/02_Variant_SNPs/core_SNP/clean.full.aln       # core alignment
results/04_CNV/cnv_depth_matrix.tsv                   # samples × regions log2 depth matrix
results/04_CNV/{sample}/cnv_per_region.tsv
results/04_Phylogeny/iqtree/iqtree.treefile
```

### Requirements

- Linux HPC with **SLURM** and **Singularity/Apptainer**
- Conda/Mamba

```bash
conda env create -f Extract_Features_Pipeline/Burkholderia_ML.yml
conda activate snakemake
```

### Configuration

The `config/` folder is git-ignored, so create it in `Extract_Features_Pipeline/` with three files:

**`config/config.yaml`**
```yaml
samples:            # filled automatically by run_pipeline.sh
  name: []
  raw_reads_r1: []
  raw_reads_r2: []
read_type: illumina
resources:
  general:
    runtime: 240    # minutes
    mem_mb: 16000
    cpus: 8
cnv_depth:
  promoter_window: 200
```

**`config/databases.yaml`**
```yaml
databases_bakta:
  bakta_path: databases/bakta_db
```

**`config/containers.yaml`** – one entry per tool (`bakta`, `seqkit`, `bwa`, `samtools`, `pilon`, `ncbi`, `fastp`, `fastqc`, `multiqc`, `bcftools`, `python_env`, `snippy`, `snpEff`, `picard`, `pindel`, `bedtools`, `mosdepth`, `freebayes`, `ariba`, `gubbins`, `iqtree`):
```yaml
snippy:
  image: docker://staphb/snippy:<version>
  output: snippy.sif
# …
```

Place the reference genome at `databases/genomes/refgenome/ref_genome.fasta` (it is annotated by Bakta during the run).

### Running

```bash
cd Extract_Features_Pipeline

# 1. Get containers
bash pull_containers.sh

# 2. Build the sample sheet (expects <sample>_R1.fastq / <sample>_R2.fastq)
#    edit reads_dir and output inside the script first
bash create_samplesheet.sh

# 3. Launch on SLURM
bash run_pipeline.sh -s samples.tsv
```

> `run_pipeline.sh` contains cluster-specific settings (`--exclude` node, `-B` bind paths). Adapt the `snakemake` call to your own HPC before running.

---

## 2. Machine-learning notebooks

Both analyses follow the same structure and should be run in numerical order.

| Step | `snp_ML/` | `cnv_ML/` |
|---|---|---|
| Exploration & cleaning | `01_data_exploration_cleaning` | `01_data_exploration_cleaning` |
| ST-aware train/test split | `02_splitting_data` | `02_splitting_data` |
| Feature selection (variance filter, SelectKBest, outlier removal) | `03_feature_selection` | `03_feature_selection` |
| PCA & K-means | `03_PCA` | `04_PCA` |
| Logistic regression | `04_logreg` / `06_logreg` | `05_logreg` |
| KNN | `05_knn_classification` | `06_knn_classification` |
| Decision tree | `05_…` / `07_decision_trees_classification` | `07_decision_trees_classification` |
| SVM | `04_…` / `06_svm_classification` | `08_svm_classification` |
| XGBoost | `07_…` / `08_xgboost_classification` | `09_xgboost_classification` |
| Nested cross-validation | `08_…` / `09_nested_cross_validation` | `10_nested_cross_validation` |
| Publication figures | `09_figures` | – |

Each model notebook trains a default baseline, tunes hyperparameters, and reports accuracy, ROC/AUC, confusion matrices and learning curves. Interpretable models also export candidate biomarkers (`cnv_ML/Results/biomarkers_*.csv`) and iTOL annotations to map predictions onto the phylogeny (`cnv_ML/Figures/*_tree.pdf`).

### Python dependencies

```bash
pip install pandas numpy scipy scikit-learn xgboost matplotlib seaborn matplotlib-venn prince jupyter
```

### Results (CNV features)

Nested cross-validation with 10 ST-separated outer test folds (`cnv_ML/nested_cv_results.csv`):

| Model | Accuracy | AUC |
|---|---|---|
| Logistic regression | 0.937 ± 0.077 | 0.971 ± 0.061 |
| XGBoost | 0.936 ± 0.082 | 0.975 ± 0.053 |
| KNN | 0.900 ± 0.073 | 0.979 ± 0.041 |

Per-model hold-out results and hyperparameters are in `cnv_ML/Results/model_results.csv`.

---

## 3. Benchmarking

`Benchmarking/AbR_reports/` contains per-isolate outputs (`*.AbR_output.final.txt` and `*_report.html`) from a rule-based, curated-mutation AMR predictor for *B. pseudomallei*, used as a reference point for comparing the ML predictions.

---

## Author

Ophélie Rutschmann
