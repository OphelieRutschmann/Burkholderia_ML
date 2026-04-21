#!/bin/bash
# Load configuration
CONFIG_FILE="config/containers.yaml"

# Snakemake workflow with Singularity/BUSCO
# Script to handle singularity image pulling and workflow execution

# Set Singularity directories to avoid space issues
export SINGULARITY_TMPDIR=$HOME/singularity_tmp
export SINGULARITY_CACHEDIR=$HOME/singularity_cache

CONT_DIR="workflow/containers"

# Create directories if they don't exist
mkdir -p $SINGULARITY_TMPDIR $SINGULARITY_CACHEDIR

#create directory where the images will be stored:
mkdir -p $CONT_DIR

echo "Singularity temporary directory: $SINGULARITY_TMPDIR"
echo "Singularity cache directory: $SINGULARITY_CACHEDIR"
echo ""

declare -A IMAGE OUTPUT

# Parse YAML file with Python
while IFS="|" read -r key image output; do
    IMAGE[$key]="$image"
    OUTPUT[$key]="$output"
done < <(
python3 - <<EOF
import yaml
import sys

with open("$CONFIG_FILE") as f:
    data = yaml.safe_load(f)

for name, props in data.items():
    print(f"{name}|{props['image']}|{props['output']}")
EOF
)

echo "Pulling bakta image"
singularity pull "$CONT_DIR/${OUTPUT[bakta]}" "${IMAGE[bakta]}"

echo "Pulling seqkit image"
singularity pull "$CONT_DIR/${OUTPUT[seqkit]}" "${IMAGE[seqkit]}"

echo "Pulling bwa image"
singularity pull "$CONT_DIR/${OUTPUT[bwa]}" "${IMAGE[bwa]}"

echo "Pulling samtools image"
singularity pull "$CONT_DIR/${OUTPUT[samtools]}" "${IMAGE[samtools]}"

echo "Pulling pilon image"
singularity pull "$CONT_DIR/${OUTPUT[pilon]}" "${IMAGE[pilon]}"

echo "Pulling ncbi-genome-download image"
singularity pull "$CONT_DIR/${OUTPUT[ncbi]}" "${IMAGE[ncbi]}"

echo "Pulling fastp image"
singularity pull "$CONT_DIR/${OUTPUT[fastp]}" "${IMAGE[fastp]}"

echo "Pulling fastqc image"
singularity pull "$CONT_DIR/${OUTPUT[fastqc]}" "${IMAGE[fastqc]}"

echo "Pulling multiqc image"
singularity pull "$CONT_DIR/${OUTPUT[multiqc]}" "${IMAGE[multiqc]}"

echo "Pulling bcftools image"
singularity pull "$CONT_DIR/${OUTPUT[bcftools]}" "${IMAGE[bcftools]}"

echo "Building python_env image"
singularity build "$CONT_DIR/${OUTPUT[python_env]}" "${IMAGE[python_env]}"

echo "Pulling snippy image"
singularity pull "$CONT_DIR/${OUTPUT[snippy]}" "${IMAGE[snippy]}"

echo "Pulling snpEff image"
singularity pull "$CONT_DIR/${OUTPUT[snpEff]}" "${IMAGE[snpEff]}"