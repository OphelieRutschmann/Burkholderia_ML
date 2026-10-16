#!/usr/bin/env bash
# Run the CNV branch of the pipeline on a local machine (no SLURM).
#
# Builds only what results/04_CNV/cnv_depth_matrix.tsv needs:
#
## The SNP, structural-variant and phylogeny branches are NOT built: this script
# names the CNV matrix as an explicit target, so Snakemake ignores `rule all`.
#
#   bash run_snakemake_local.sh -s samples.tsv [-j N] [-m MB] [-n] [-u]
#     -s  sample sheet (sample_name, raw_reads_r1, raw_reads_r2)   [required]
#     -j  cores to use                        (default: all but one)
#     -m  memory ceiling in MB                (default: 80% of RAM)
#     -n  dry run - print the DAG, build nothing
#     -u  unlock a working directory left locked by a killed run
#
# Run it from Extract_Features_Pipeline/, the same place as run_pipeline.sh.
set -uo pipefail

SHEET=""; CORES=""; MEM=""; DRY=0; UNLOCK=0
while getopts "s:j:m:nuh" opt; do
  case $opt in
    s) SHEET="$OPTARG" ;;
    j) CORES="$OPTARG" ;;
    m) MEM="$OPTARG" ;;
    n) DRY=1 ;;
    u) UNLOCK=1 ;;
    h) sed -n '2,20p' "$0"; exit 0 ;;
    *) exit 1 ;;
  esac
done
[[ -z "$SHEET" ]] && { echo "ERROR: -s sample sheet is required" >&2; exit 1; }
[[ -f "$SHEET" ]] || { echo "ERROR: no such sample sheet: $SHEET" >&2; exit 1; }

# ---------------------------------------------------------------- preflight
fail=0
note() { echo "  [ok]   $*"; }
bad()  { echo "  [FAIL] $*" >&2; fail=1; }

echo "==> checking the local environment"

command -v snakemake >/dev/null && note "snakemake $(snakemake --version 2>/dev/null)" \
  || bad "snakemake not on PATH - conda activate your pipeline env first"

CONTAINER_FLAG=""
if command -v apptainer >/dev/null; then
  CONTAINER_FLAG="--use-apptainer"; note "apptainer $(apptainer --version 2>/dev/null | awk '{print $NF}')"
elif command -v singularity >/dev/null; then
  CONTAINER_FLAG="--use-singularity"; note "singularity $(singularity --version 2>/dev/null | awk '{print $NF}')"
else
  bad "neither apptainer nor singularity found - the rules all declare containers,
         so there is no conda fallback. Install Apptainer, then: bash pull_containers.sh"
fi

[[ -f config/config.yaml    ]] && note "config/config.yaml"    || bad "config/config.yaml missing (config/ is gitignored - see the README)"
[[ -f config/databases.yaml ]] && note "config/databases.yaml" || bad "config/databases.yaml missing"
[[ -f config/containers.yaml ]] && note "config/containers.yaml" || bad "config/containers.yaml missing"

REF="databases/genomes/refgenome/ref_genome.fasta"
[[ -f "$REF" ]] && note "reference genome" || bad "$REF missing - put the reference there first"

# The CNV branch needs these images. fastqc is not needed: the CNV target never
# asks for its output, only fastp's JSON reaches multiqc.
NEED=(fastp multiqc bwa samtools picard bedtools mosdepth python)
missing=()
for c in "${NEED[@]}"; do [[ -f "workflow/containers/${c}.sif" ]] || missing+=("${c}.sif"); done
if [[ ${#missing[@]} -gt 0 ]]; then
  bad "missing containers: ${missing[*]}
         run: bash pull_containers.sh
         (note: the CNV rules ask for python.sif, while containers/python_env.def
          builds python_env.sif - check what containers.yaml writes it as)"
else
  note "all ${#NEED[@]} CNV containers present"
fi

# Bakta only matters if the annotation has not been produced yet. If you already
# have ref_genome.gbff from the HPC run, copy it in and skip the database entirely.
GBFF="databases/genomes/refgenome/ref_genome.gbff"
if [[ -f "$GBFF" ]]; then
  note "ref_genome.gbff present - Bakta will be skipped"
elif [[ -f "workflow/containers/bakta.sif" ]]; then
  note "ref_genome.gbff absent - Bakta will run (downloads a ~1.3 GB light DB, annotation takes 10-30 min)"
else
  bad "ref_genome.gbff absent and no bakta.sif - either copy the .gbff from the HPC
         into databases/genomes/refgenome/, or pull the Bakta container"
fi

# Known bug: the rule calls gbff_to_gff3.py, the file on disk is gbff_gff3.py.
if [[ ! -e workflow/scripts/gbff_to_gff3.py ]]; then
  if [[ -f workflow/scripts/gbff_gff3.py ]]; then
    ln -s gbff_gff3.py workflow/scripts/gbff_to_gff3.py
    echo "  [FIXED] variant_calling_cnvs.smk calls workflow/scripts/gbff_to_gff3.py,"
    echo "          but the file is named gbff_gff3.py. Created a symlink so the run"
    echo "          proceeds. The real fix is to correct the name in the .smk file."
  else
    bad "neither gbff_to_gff3.py nor gbff_gff3.py found in workflow/scripts/"
  fi
fi

[[ $fail -eq 1 ]] && { echo; echo "preflight failed - fix the above and rerun" >&2; exit 1; }

# ---------------------------------------------------------------- resources
if [[ -z "$CORES" ]]; then
  total=$( (command -v nproc >/dev/null && nproc) || sysctl -n hw.ncpu 2>/dev/null || echo 4 )
  CORES=$(( total > 1 ? total - 1 : 1 ))
fi
if [[ -z "$MEM" ]]; then
  if [[ -r /proc/meminfo ]]; then
    MEM=$(awk '/MemTotal/ {printf "%d", $2/1024*0.8}' /proc/meminfo)
  elif command -v sysctl >/dev/null; then
    MEM=$(( $(sysctl -n hw.memsize) / 1048576 * 8 / 10 ))
  else
    MEM=8000
  fi
fi
echo "==> running with $CORES cores, ${MEM} MB memory ceiling"

# ---------------------------------------------------------------- config
cp config/config.yaml "config/config.yaml.bak.$(date +%Y%m%d_%H%M%S)"
python3 - "$SHEET" <<'PY' || exit 1
import csv, sys, yaml
sheet = sys.argv[1]
with open(sheet) as f:
    rows = [r for r in csv.DictReader(f, delimiter="\t")]
if not rows:
    sys.exit("ERROR: sample sheet is empty")
missing = [r["sample_name"] for r in rows
           if not r.get("raw_reads_r1") or not r.get("raw_reads_r2")]
if missing:
    sys.exit("ERROR: these samples have no read paths: " + ", ".join(missing))
with open("config/config.yaml") as f:
    cfg = yaml.safe_load(f)
cfg.setdefault("samples", {})
cfg["samples"]["name"] = [r["sample_name"] for r in rows]
cfg["samples"]["raw_reads_r1"] = [r["raw_reads_r1"] for r in rows]
cfg["samples"]["raw_reads_r2"] = [r["raw_reads_r2"] for r in rows]
with open("config/config.yaml", "w") as f:
    yaml.dump(cfg, f, default_flow_style=False, sort_keys=False)
print(f"==> {len(rows)} samples written to config/config.yaml")
PY

# Every directory holding reads has to be visible inside the containers.
BINDS="-B $(pwd)"
while IFS=$'\t' read -r _ r1 _; do
  d=$(dirname "$r1"); [[ -d "$d" ]] && BINDS="$BINDS -B $d"
done < <(tail -n +2 "$SHEET" | awk -F'\t' '{print $1"\t"$2"\t"$3}')
BINDS=$(tr ' ' '\n' <<<"$BINDS" | awk '!seen[$0]++' | tr '\n' ' ')

export TMPDIR="$(pwd)/tmp_scratch"; mkdir -p "$TMPDIR"
export APPTAINER_TMPDIR="$(pwd)/singularity_tmp"
export SINGULARITY_TMPDIR="$APPTAINER_TMPDIR"; mkdir -p "$APPTAINER_TMPDIR"
BINDS="$BINDS -B ${TMPDIR}:/tmp"
trap 'rm -rf "$APPTAINER_TMPDIR" "$TMPDIR"' EXIT

if [[ $UNLOCK -eq 1 ]]; then
  echo "==> unlocking"; snakemake --unlock --cores 1; exit $?
fi

# ---------------------------------------------------------------- run
# Naming the CNV matrix explicitly means `rule all` (SNPs, SVs, phylogeny) is
# never evaluated. The per-sample cnv_per_region.tsv files are inputs to it.
TARGET="results/04_CNV/cnv_depth_matrix.tsv"

echo "==> target: $TARGET"
[[ $DRY -eq 1 ]] && echo "==> dry run"

snakemake "$TARGET" \
    --cores "$CORES" \
    --resources mem_mb="$MEM" \
    $CONTAINER_FLAG \
    --singularity-args "$BINDS" \
    --latency-wait 30 \
    --rerun-incomplete \
    --printshellcmds \
    --keep-going \
    $( [[ $DRY -eq 1 ]] && echo "--dry-run --quiet" )
rc=$?

echo
if [[ $rc -eq 0 && $DRY -eq 0 ]]; then
  echo "Pipeline finished."
  if [[ -f "$TARGET" ]]; then
    echo "  matrix : $TARGET"
    echo "  samples: $(( $(wc -l < "$TARGET") - 1 ))"
    echo
    echo "Samples that failed QC are listed out of results/00_QC/passed_samples.csv;"
    echo "compare it against your sample sheet to see what was dropped."
  fi
elif [[ $rc -ne 0 ]]; then
  echo "Pipeline failed (exit $rc). --keep-going was set, so unrelated samples may" >&2
  echo "still have completed; scroll up for the first real error." >&2
fi
exit $rc