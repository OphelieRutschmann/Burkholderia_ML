#!/bin/bash
# Script to run the pipeline
# Usage: ./run_pipeline.sh -s samples.tsv

# Default values
SAMPLE_SHEET=""

# Parse options
while getopts ":s:" opt; do
  case $opt in
    s)
      SAMPLE_SHEET="$OPTARG"
      ;;
    \?)
      echo "Invalid option: -$OPTARG" >&2
      echo "Usage: $0 -s samples.tsv"
      exit 1
      ;;
    :)
      echo "Option -$OPTARG requires an argument." >&2
      echo "Usage: $0 -s samples.tsv"
      exit 1
      ;;
  esac
done

# Check if sample sheet provided
if [ -z "$SAMPLE_SHEET" ]; then
    echo "ERROR: Sample sheet must be provided with -s"
    echo "Usage: $0 -s samples.tsv"
    exit 1
fi

echo "Running pipeline with:"
echo "  Sample sheet: $SAMPLE_SHEET"
echo ""

# Update config with samples from sample sheet
python3 << EOF
import yaml
import csv
import sys

# Read sample sheet
samples = []
r1_paths = []
r2_paths = []

with open('$SAMPLE_SHEET', 'r') as f:
    reader = csv.DictReader(f, delimiter='\t')
    rows = list(reader)
    
    if not rows:
        print("ERROR: Sample sheet is empty!")
        sys.exit(1)
    
    for row in rows:
        samples.append(row["sample_name"])
        r1_paths.append(row.get("raw_reads_r1", ""))
        r2_paths.append(row.get("raw_reads_r2", ""))
       
# Update config
with open('config/config.yaml', 'r') as f:
    config = yaml.safe_load(f)

config['samples']['name'] = samples
config['samples']['raw_reads_r1'] = r1_paths
config['samples']['raw_reads_r2'] = r2_paths

# Write updated config
with open('config/config.yaml', 'w') as f:
    yaml.dump(config, f, default_flow_style=False, sort_keys=False)

print(f"Loaded {len(samples)} samples: {', '.join(samples)}")
EOF

# Make Singularity temp dir with absolute path
export SINGULARITY_TMPDIR="$(pwd)/singularity_tmp"
mkdir -p "$SINGULARITY_TMPDIR"
echo "Using Singularity tmpdir: $SINGULARITY_TMPDIR"
df -h "$SINGULARITY_TMPDIR" | tail -1
echo ""

# Cleanup temp dir at end of run
cleanup_tmpdir() {
    echo ""
    echo "Cleaning up Singularity tmpdir: $SINGULARITY_TMPDIR"
    rm -rf "$SINGULARITY_TMPDIR"
    echo "Cleanup complete"
}
# Register cleanup to run on exit (success or failure)
trap cleanup_tmpdir EXIT

# Run snakemake
echo ""
echo "Starting Snakemake..."
snakemake \
    --cluster "sbatch --cpus-per-task={resources.cpus_per_task} --mem={resources.mem_mb} --time={resources.runtime} --exclude=hpc-srvbio-01" \
    --jobs 10 \
    --resources mem_mb=400000 \
    --use-singularity \
    --use-conda --latency-wait=60 \
    --conda-frontend conda \
    --rerun-incomplete \
    --singularity-args "-B /mnt/nfs -B /home/ruop/" \

if [ $? -eq 0 ]; then
    echo ""
    echo "Pipeline completed successfully!"
else
    echo ""
    echo "Pipeline failed!"
    exit 1
fi