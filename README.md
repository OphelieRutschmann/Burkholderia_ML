# Burkholderia_ML
Predicting co-trimoxazole (trimethoprim–sulfamethoxazole, **SXT**) resistance in *Burkholderia pseudomallei* from whole-genome sequencing data using machine learning.

## Author
Ophélie Hofer-Rutschmann

## Introduction
*Burkholderia pseudomallei* is a Gram-negative bacterium and the causative agent of melioidosis, a disease with a mortality rate exceeding 50%. Treatment is complicated by the bacterium’s many intrinsic and acquired antibiotic resistances, and determining the resistance profile of a clinical strain is a major bottleneck. Current laboratory-based susceptibility tests are slow and labor-intensive, and sequencing combined with *in silico* resitance prediction is largely ineffective for *B. pseudomallei*. Indeed, unlike most bacteria, *B. pseudomallei* develops resistances through gene copy number variations (CNVs) and single-nucleotide polymorphisms (SNPs) rather than through the acquisition of new resistance genes that could be detected using existing *in silico* resistance prediction tools.

This project takes advantage of ~300 publically available *B. pseudomallei* genomes with corresponding resistance profiles to assess the feasability of machine-learning-based approaches for resistance prediction. To this end, we developed a Snakemake Pipeline to extract genomic features from the raw sequencing data that can be used as input for a machine learning model. Three kinds of feature are extracted across both the core and the accessory genome: copy number variations (CNVs), single nucleotide polymorphisms (SNPs) and structural variants (SVs). SVs complement the other two: a CNV is seen as a change in read depth over a region, whereas an SV is called from the reads themselves and gives the breakpoints of a deletion, a tandem duplication or an inversion.

## Contents
The repository has two parts:

1. **`Extract_Features_Pipeline/`**: This is a Snakemake Pipeline extracts genomic features from raw paired-end Illumina reads. The Pipeline extracts core-genome variants (SNPs, small deletions / insertions), copy number variations in gene and promoter regions and structural variants (deletions, tandem duplications and inversions), and also performs a core-SNP phylogeny.
2. **`snp_ML/`, `cnv_ML/` and `sv_ML/`**: These files contain Jupyter notebooks that clean the feature matrices, select features and train/compare classifiers (logistic regression, KNN, SVM, decision tree, XGBoost). The three folders follow the same structure, so the three feature types can be compared directly. The final model performance is estimated by lineage–aware nested cross-validation.

---

## Dataset
For this project, publicly available raw sequencing data from the following publications was used:

- Hinwan Y, Chareonsudjai P, Reungsang P, Kraiklang R, Chetchotisakd P, Chareonsudjai S, Sirichoat A, Nithichanon A, Wonglakorn L, Sermswan RW, Blair D, Faksri K. Analysis of fine-scale phylogeny of Burkholderia pseudomallei in relation to regional geography and drug susceptibility in Thailand. Sci Rep. 2024 Aug 28;14(1):19961. doi: 10.1038/s41598-024-70558-5.
- Viberg LT, Sarovich DS, Kidd TJ, Geake JB, Bell SC, Currie BJ, Price EP. Within-Host Evolution of Burkholderia pseudomallei during Chronic Infection of Seven Australasian Cystic Fibrosis Patients. mBio. 2017 Apr 11;8(2):e00356-17. doi: 10.1128/mBio.00356-17.
- Madden DE, Webb JR, Steinig EJ, Currie BJ, Price EP, Sarovich DS. Taking the next-gen step: Comprehensive antimicrobial resistance detection from Burkholderia pseudomallei. EBioMedicine. 2021 Jan;63:103152. doi: 10.1016/j.ebiom.2020.103152. Epub 2020 Dec 4.
- Bugrysheva JV, Sue D, Hakovirta J, Loparev VN, Knipe K, Sammons SA, Ranganathan-Ganakammal S, Changayil S, Srinivasamoorthy G, Weil MR, Tatusov RL, Gee JE, Elrod MG, Hoffmaster AR, Weigel LM. Finished Annotated Genome Sequence of Burkholderia pseudomallei Strain Bp1651, a Multidrug-Resistant Clinical Isolate. Genome Announc. 2015 Dec 3;3(6):e01427-15. doi: 10.1128/genomeA.01427-15.

Alltogether, the data consists of **282 genomes** with corresponding antibiotic resistance profiles for 5 antibiotics (Co-Trimexazole (SXT), Ceftazidime (CAZ), Amoxicillin–clavulanic acid (AMC), Meropenem (MEM) and Imipenem (IPM)). For this project, we will focus on predicting resistance to **Co-Trimexazole** only. 

The raw sequencing reads was downloaded from SRA and relevant features from the data were extracted using a custom Snakemake pipeline. 

## Repository structure

```
Burkholderia_ML/
├── Extract_Features_Pipeline/      # Snakemake feature-extraction workflow
│   ├── Snakefile
│   ├── run_pipeline.sh             # launches the pipeline on SLURM
│   ├── create_samplesheet.sh       # builds samples.tsv from the raw reads folder
│   ├── pull_containers.sh          # pulls/builds Singularity images
│   ├── Burkholderia_ML.yml         # conda env (Snakemake 7.32)
│   └── workflow/
│       ├── modules/                # QC, databases, SNPs, SVs, CNVs, phylogeny
│       ├── rules/common.smk        # shared helper functions
│       ├── scripts/                # Python helpers (matrix conversion, normalisation, etc.)
│       └── containers/python_env.def
├── snp_ML/                         # ML on core-genome SNP features
│   ├── Data/                       
│   ├── Figures/
│   ├── Results/ 
│   └── 01_… – 10_….ipynb           # Analysis notebooks
├── cnv_ML/                         # ML on copy-number features
│   ├── Data/                       
│   ├── Figures/                    
│   ├── Results/                    
│   ├── nested_cv_results.csv
│   ├── 01_… – 10_….ipynb           # Analysis notebooks
│   ├── 11_binarised_cnv.ipynb      # the same models on presence/absence CNV calls
│   └── biomarker_comparison.ipynb
├── sv_ML/                          # ML on structural-variant features
│   ├── Data/                       
│   ├── Figures/                    
│   ├── Results/                    
│   ├── 01_… – 10_….ipynb           # Analysis notebooks
│   └── biomarker_comparison.ipynb
├── genome_metadata*.csv
└── itol_SXT
```

---

## 1. Feature extraction pipeline

### Description
For this project,a Snakemake pipeline was developed to extract copy number variations (CNVs), genomic variants (mainyl SNPs, but also small insertions and deletions) and structural variants (SVs) from the genome of each strain. **CNVs** are changes in the genome that result in the *deletion* or *duplication* of entire DNA regions (genes, promoters or non-coding regions). **SNPs** are genetic variations where a single nucleotide differs at a specific position in the DNA sequence between strains. **SVs** are rearrangements of at least 50 bp reported with their breakpoints: deletions, tandem duplications and inversions. These genetic features were chosen for the project as CNVs, small genetic variants and genome rearrangements have all been linked to the acquisition of antibiotic resistances in *B. pseudomallei*. CNVs and SVs describe overlapping events from two different angles: a CNV is a change in read depth over an annotated region, while an SV is called from the read pairs and split reads and gives the position of the breakpoints, which also makes inversions visible (they change no copy number).

The Snakemake pipeline automatizes the extraction of these features from the genomes and contains the following modules:

- **databases.smk**: this module setups the Bakta database necessary for the pipeline.
- **QC_illumina.smk**: this module performs quality control, read preprocessing and filtering. It runs **FastQC** on the raw sequencing data, and performs trimming and preliminary quality filtering with **fastp**. QC metrics are aggregated with **MultiQC**, and the trimmed reads are aligned to a reference genome using **BWA-MEM**. Coverage and mapping statistics are calculated with **samtools**, and a filtering step is performed to remove samples that fail any of the following criteria post trimming: a Q30 rate >90%, a GC content within ±3% of the expected 67%, at least 75% of reads surviving fastp trimming/filtering, and a mean genome coverage above 1x. The list of samples that pass this filtering step are written in *passed_samples.csv* which is then used by the downstream modules.
- **variant_calling_snps.smk**: this module extract core-SNPs (SNPs in the core genome) for each sample. First, variants are called in comparison to a reference genome using **Snippy**. The reference genome is then annotated with **Bakta** and a **SnpEff** database is built using this annotation. A core-genome SNP alignment is constructed using all samples (with **snippy-core**) to identify SNPs present in the core genome (the genome that is common to most samples). This core-genome SNP alignment is then converted into a presence/absence SNP feature matrix, and the SNPs are annotated with SnpEff to classify them by their predicted impact (LOW, MODERATE, HIGH, MODIFIER) and effect (disruptive_inframe_deletion, missense_variant, upstream_gene_variant, etc.) on the gene they fall in.
- **variant_calling_cnvs.smk**: this module calls copy-number variations (CNVs) for each annotated gene / region in the genome using read depth. First, the annotated reference genome (GBFF) generated in the previous module is converted to GFF3 (necessary for the analysis). Then, a BED file containing the genes with their upstream promoter regions is built. Per-region sequencing depth is then computed using **mosdepth** using the deduplicated BAM files generated with **samtools** in the **QC_illumina.smk** module. The copy number (the depth) is normalized for each sample for each chromosome by dividing each region's raw coverage by that chromosome's median depth. The output of this step is a cnv_ratio where 1 suggests typical coverage (1 copy of the region is present), and values above or below 1 suggest a duplication or deletion. This ratio is then log2-transformed (0 indicates no change, while positive/negative values indicate duplications/deletions. Absent regions are set to NaN), and the results are merged into a CNV depth matrix. Contrary to the variant_calling_snps.smk module, CNVs are computed in both the core genome (common to most strains) and the accessory genome. The output of this module is the file *cnv_depth_matrix.tsv* which contains the log2 transformed depth for each annotated region for each strain.
- **variant_calling_svs.smk**: this module calls structural variants (SVs) for each sample. PCR duplicates are first marked with **Picard**, then **Pindel** is run on the deduplicated alignments to call deletions, tandem duplications and inversions. Long insertions are left out, as Pindel reports them without a size. Per-sample calls are converted to VCF with pindel2vcf, concatenated and sorted with **bcftools**, and then merged across samples with **SURVIVOR**, which treats two calls as the same event when their breakpoints lie within `max_dist` bp of one another. This tolerance matters: Pindel breakpoints shift by a few bases between strains, so merging on exact positions would split one real event into many features. Only SVs of at least `min_size` bp are kept.

  An SV is not called in a sample for two very different reasons: either the region is intact, or the region is missing from that strain altogether. To tell these apart, the flanks of every merged SV are checked for read coverage in every sample with **mosdepth**. A flank counts as covered when at least `min_frac` of its bases reach `min_depth`. Each cell of the matrix is then **1** (SV called), **0** (not called but both flanks covered, so the SV is genuinely absent) or **-1** (not called and a flank is missing, so the locus cannot be scored in that strain). An SV is kept as a *core* SV when its locus is covered in at least `min_covered_samples` of the samples; the remaining -1 values are then set to 0. SVs called in all samples or in none carry no information and are dropped. The outputs of this module are *svs_core.tsv* and *svs_core_meta.tsv*.

- **phylogeny.smk**: this module uses the core-SNPs to generate a phylogenetic tree using **iqtree2**. This tree shows us how the strains are related to each other and more importantly, how the resistance profile is distributed according to phylogeny. 


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
| SV calling | `variant_calling_svs.smk` | Picard MarkDuplicates, Pindel (deletions, tandem duplications, inversions), bcftools, SURVIVOR (`max_dist` 100 bp, `min_size` 50 bp) |
| SV locus coverage & matrix | `variant_calling_svs.smk`, `convert_sv_to_matrix.py` | mosdepth over ±500 bp flanks; SV kept when its locus is covered in ≥ 95 % of samples |
| Phylogeny | `phylogeny.smk` | IQ-TREE 2 (GTR+G) on the core alignment |

All tools run in Singularity/Apptainer containers stored in `workflow/containers/`.

### Main outputs

The pipeline outputs five .tsv files that will be used to train machine learning models: 

- **cnv_depth_matrix.tsv**: this files contains the CNVs for all identified genomic regions for each sample. As described above, the CNVs are computed as *log2(normalized depth for the region of interest)*. 0 indicates no change in that specific region (the region is present in one copy). Negative or positive values indicate deletions or duplications, respectively. Absent regions are set to NaN. The feature matrix dimensions are n_samples x n_genomic_region.
- **snps.tsv**: this file contains the SNPs identified for each strain. This feature matrix is one-hot encoded: each possible SNP is a feature that can be either absent (0) or present (1) for a specific strain. Example: the feature snp__NZ_AP028077.1_128_G_C, corresponding to a G to C mutation in chromosome NZ_AP028077.1 at position 128, can either be present (1), signifying that the strain has that mutation, or absent (0), signifying that the strain does not have that mutation. As mentioned above, only SNPs in the core genome (regions present in most of the strains) are included in the analysis. The feature matrix dimensions are n_samples x n_identified_core_SNPs. 
- **svs_core.tsv**: this file contains the presence/absence of each core structural variant for each sample, coded 1 (the SV was called in that strain) or 0 (it was not). Feature names follow *sv__{chromosome}_{start}_{end}_{type}*, for example sv__NZ_CP009156.1_1254_1350_DEL for a 96 bp deletion on chromosome NZ_CP009156.1. The feature matrix dimensions are n_samples x n_core_SVs. Only SVs whose locus could be scored in nearly all samples are included, so a 0 means the region was examined and found intact rather than missing from the assembly (see the description of variant_calling_svs.smk above).
- **svs_core_meta.tsv**: this file contains information about each SV in svs_core.tsv: its type (DEL, DUP or INV), its length in bp, the number of samples carrying it and the fraction of samples in which its locus was covered. It can be used to filter SVs by type or size for downstream analysis, in the same way snps_meta.tsv is used for the SNPs.
- **snps_meta.tsv**: this file was generated using SnpEff and contains information about all of the SNPs present in the snps.tsv file. For each SNP, the corresponding gene, predicted impact (LOW, MODERATE, HIGH, MODIFIER) and predicted effect (disruptive_inframe_deletion, missense_variant, upstream_gene_variant, conservative_inframe_deletion, etc.) are reported. This file can be used to filter out SNPs for downstream analysis.

In addition to these files, a **genome_metadata.csv** file containing information about the strains (genome_name, source, species, resistance_profile, sra_accession, Sequence Type, location) will also be used for downstream analysis. This metadata file was generated manually by summarizing the information contained in the publications cited above.

### Requirements

- Linux HPC with **SLURM** and **Singularity/Apptainer**
- Conda/Mamba

```bash
conda env create -f Extract_Features_Pipeline/Burkholderia_ML.yml
conda activate snakemake
```

### Configuration

`Extract_Features_Pipeline/config/` holds three files. They are included in the repository, so a fresh clone can run the pipeline; the paths inside them are cluster-specific and should be adapted. The `.gitignore` entry for `config/*.yaml` only stops *new* YAML files from being added by accident.

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
sv:                            # all optional, the values below are the defaults
  max_dist: 100                # max breakpoint distance (bp) to merge calls across samples
  min_size: 50                 # minimum SV length (bp)
  flank: 500                   # flank size (bp) checked on each side of an SV
  min_depth: 10                # depth a flank base needs to count as covered
  min_frac: 0.9                # fraction of flank bases that must reach min_depth
  min_covered_samples: 0.95    # fraction of samples in which an SV locus must be covered
```

> `min_covered_samples` is the setting to watch. Each strain covers roughly 90 % of the SV loci, but a *different* 90 % per strain, so requiring every locus to be covered in all samples keeps almost nothing. The run log prints how many SVs would survive at 1.0, 0.99, 0.95, 0.9 and 0.8, together with the samples whose SV loci are least covered, so the threshold can be chosen from the data.

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

### Running the pipeline

```bash
cd Extract_Features_Pipeline

# 1. Get containers
bash pull_containers.sh

# 2. Build the sample sheet (expects <sample>_R1.fastq / <sample>_R2.fastq)
#    edit reads_dir and output inside the script first
bash create_samplesheet.sh

# 3. Run on SLURM
bash run_pipeline.sh -s samples.tsv
```

> `run_pipeline.sh` contains cluster-specific settings (`--exclude` node, `-B` bind paths). Adapt the `snakemake` call to your own HPC before running.

---

## 2. Machine-learning notebooks

The SNP, CNV and SV analyses follow the same structure and should be run in numerical order.

| Step | Notebook |
|---|---|
| Exploration & cleaning | `01_data_exploration_cleaning` |
| ST-aware train/test split | `02_splitting_data` |
| Feature selection (variance filter, SelectKBest, outlier removal) | `03_feature_selection` |
| PCA & K-means | `04_PCA` |
| Logistic regression | `05_logreg` |
| KNN | `06_knn_classification` |
| Decision tree | `07_decision_trees` |
| SVM | `08_svm` |
| XGBoost | `09_xgboost` |
| Nested cross-validation | `10_nested_cross_validation` |
| Comparison of the biomarkers found by the different models | `biomarker_comparison` (`cnv_ML`, `sv_ML`) |
| The same models on presence/absence CNV calls | `11_binarised_cnv` (`cnv_ML` only) |

In each of the model-specific notebook, a basic model is trained, then hyperparameters are tuned and the accuracy, ROC/AUC, confusion matrices and learning curves are calculated. Candidate biomarkers are also extracted from the best models and stored in (`/Results/biomarkers_*.csv`)

The feature matrices of the three folders differ in what a feature is, so the notebooks differ in a few places. SNPs and SVs are binary, so they are filtered by variance and scored with a chi2 test, and features sharing a presence/absence pattern are collapsed into one. CNVs are continuous log2 depths, so they are filtered by standard deviation, scored with an ANOVA F-test and grouped with `FeatureAgglomeration`. For the SVs, `03_feature_selection` keeps the 200 best features by chi2 score rather than only the significant ones, since almost no single SV reaches significance after correction for multiple testing.

`10_nested_cross_validation` repeats the feature selection inside every outer fold, so no fold takes part in choosing the features it is tested on, and each fold holds whole sequence types. `11_binarised_cnv` in `cnv_ML` exists because depth at `dnaA`, the gene at the origin of replication, discriminates resistant from susceptible strains almost as well as the whole CNV matrix: the continuous depths partly measure how fast the cells were replicating when the DNA was extracted. That notebook repeats the nested cross-validation on presence/absence calls, which keeps the region-specific signal and drops most of that effect; notebook 10 is left unchanged so the two can be compared.

### Python dependencies

```bash
pip install pandas numpy scipy scikit-learn xgboost matplotlib seaborn matplotlib-venn prince jupyter
```
