#!/usr/bin/env bash
# Download the Hinwan et al. 2024 genomes (BioProject PRJNA1051349) that are not
# yet in the local feature set.
#
# These are the isolates whose SXT phenotype was "ND" in the paper's own
# supplementary table but for which the 2025 Infection, Genetics and Evolution
# paper publishes a co-trimoxazole MIC. They are therefore newly usable as
# labelled data, and none of them has been through the pipeline before.
#
# The accession list is read from a metadata TSV (default samples_new_metadata.tsv,
# 169 rows) rather than hardcoded here, so editing that file changes what is
# fetched. Reads come from the ENA mirror over HTTPS: nothing to install,
# resumable, and md5 checksums are published alongside the files.
#
#   bash download_hinwan_new.sh [-m TSV] [-o DIR] [-j N] [-e DIR] [-u] [-n]
#     -m  metadata TSV with an 'accession' column  (default: samples_new_metadata.tsv)
#     -o  output directory                        (default: raw_reads_new)
#     -j  parallel downloads                      (default: 3; be polite to ENA)
#     -e  extra directory of existing reads to check first (repeatable)
#     -u  also decompress to .fastq               (roughly 4x the disk)
#     -n  dry run: resolve URLs, print sizes, download nothing
#
# Writes samples_new.tsv in the pipeline's sample-sheet format. Re-running skips
# files already present whose md5 matches, so it is safe to resume or repeat.
#
# Expect ~114 GB compressed (372 Gbp across 169 runs; median ~292x of the 7.25 Mb
# genome, minimum 127x). Check free space before starting, and use -n first.
set -uo pipefail

META="samples_new_metadata.tsv"
OUTDIR="raw_reads_new"
JOBS=3
DECOMPRESS=0
DRYRUN=0
EXTRA_DIRS=()

while getopts "m:o:j:e:unh" opt; do
  case $opt in
    m) META="$OPTARG" ;;
    o) OUTDIR="$OPTARG" ;;
    j) JOBS="$OPTARG" ;;
    e) EXTRA_DIRS+=("$OPTARG") ;;
    u) DECOMPRESS=1 ;;
    n) DRYRUN=1 ;;
    h) sed -n '2,28p' "$0"; exit 0 ;;
    *) exit 1 ;;
  esac
done

command -v curl >/dev/null || { echo "ERROR: need curl" >&2; exit 1; }
[[ -f "$META" ]] || { echo "ERROR: no metadata at $META (use -m)" >&2; exit 1; }
MD5BIN=$(command -v md5sum || command -v md5 || true)
[[ -z "$MD5BIN" ]] && echo "WARNING: no md5sum/md5 found - downloads will not be verified" >&2

mkdir -p "$OUTDIR"
MANIFEST="$OUTDIR/.manifest.tsv"

# ---- read accession + label from the metadata, by column name ----
# Column order is not assumed; the header is parsed for 'accession' and, if
# present, 'genome_name' and 'resistance_profile' for the run log.
LISTING=$(awk -F'\t' '
  NR==1 {
    for (i=1; i<=NF; i++) { gsub(/\r/,"",$i); h[$i]=i }
    if (!("accession" in h)) { print "ERR:no accession column" > "/dev/stderr"; exit 1 }
    a=h["accession"]; g=("genome_name" in h)?h["genome_name"]:0
    p=("resistance_profile" in h)?h["resistance_profile"]:0
    next
  }
  $a != "" {
    gsub(/\r/,"",$a)
    label = (g ? $g : "-")
    if (p && $p != "") { sxt = ($p ~ /Resistant to[^;]*SXT/) ? "R" : "S"; label = label " " sxt }
    print $a "\t" label
  }
' "$META") || exit 1
[[ -n "$LISTING" ]] || { echo "ERROR: no accessions found in $META" >&2; exit 1; }

mapfile -t ACCESSIONS < <(cut -f1 <<<"$LISTING")
declare -A LABEL
while IFS=$'\t' read -r a l; do LABEL["$a"]="$l"; done <<<"$LISTING"
echo "==> $META: ${#ACCESSIONS[@]} accessions"

# ---- skip anything already on disk ----
have_reads() {
  local acc="$1" d pat
  for d in "$OUTDIR" ${EXTRA_DIRS+"${EXTRA_DIRS[@]}"}; do
    [[ -d "$d" ]] || continue
    for pat in "${acc}_R1.fastq.gz" "${acc}_R1.fastq" "${acc}_1.fastq.gz" "${acc}_1.fastq"; do
      [[ -s "$d/$pat" ]] && { echo "$d/$pat"; return 0; }
    done
  done
  return 1
}

TODO=(); SKIPPED=0
for acc in "${ACCESSIONS[@]}"; do
  if hit=$(have_reads "$acc"); then
    echo "  - $acc (${LABEL[$acc]:-}) already present at $hit"
    SKIPPED=$((SKIPPED+1))
  else
    TODO+=("$acc")
  fi
done
echo "==> $SKIPPED already on disk | ${#TODO[@]} to fetch"
[[ ${#TODO[@]} -eq 0 ]] && { echo "nothing to do"; exit 0; }

# ---- resolve each run at ENA ----
echo "==> resolving ${#TODO[@]} runs at ENA"
: > "$MANIFEST"
FAILED_RESOLVE=()
for acc in "${TODO[@]}"; do
  row=$(curl -fsS --retry 3 --retry-delay 2 \
    "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=${acc}&result=read_run&fields=run_accession,fastq_ftp,fastq_md5,fastq_bytes&format=tsv" \
    2>/dev/null | awk 'NR==2')
  if [[ -z "$row" ]]; then
    echo "  !! $acc  no ENA record" >&2; FAILED_RESOLVE+=("$acc"); continue
  fi
  ftp=$(cut -f2 <<<"$row"); md5=$(cut -f3 <<<"$row"); bytes=$(cut -f4 <<<"$row")
  if [[ "$ftp" != *";"* ]]; then
    echo "  !! $acc  not paired-end at ENA" >&2; FAILED_RESOLVE+=("$acc"); continue
  fi
  printf '%s\t%s\t%s\t%s\n' "$acc" "$ftp" "$md5" "$bytes" >> "$MANIFEST"
  mb=$(( $(cut -d';' -f1 <<<"$bytes")/1048576 + $(cut -d';' -f2 <<<"$bytes")/1048576 ))
  printf '  %-13s %6d MB  %s\n' "$acc" "$mb" "${LABEL[$acc]:-}"
done

n=$(wc -l < "$MANIFEST")
[[ "$n" -eq 0 ]] && { echo "ERROR: nothing resolved" >&2; exit 1; }
total=$(awk -F'\t' '{split($4,b,";"); s+=b[1]+b[2]} END{printf "%.1f", s/1073741824}' "$MANIFEST")
echo "==> $n runs resolved, ${total} GB compressed"
[[ ${#FAILED_RESOLVE[@]} -gt 0 ]] && echo "    unresolved: ${FAILED_RESOLVE[*]}" >&2
[[ $DECOMPRESS -eq 1 ]] && echo "    (-u set: expect roughly 4x that on disk)"
if [[ $DRYRUN -eq 1 ]]; then echo "dry run, stopping here"; exit 0; fi

# ---- fetch ----
fetch_one() {
  local acc="$1" ftp="$2" md5="$3" i=0 url want dest have
  for url in $(tr ';' ' ' <<<"$ftp"); do
    i=$((i+1)); want=$(cut -d';' -f$i <<<"$md5")
    dest="${OUTDIR}/${acc}_R${i}.fastq.gz"
    if [[ -s "$dest" && -n "$MD5BIN" ]]; then
      have=$($MD5BIN "$dest" | awk '{print $1}')
      [[ "$have" == "$want" ]] && { echo "  = ${acc}_R${i} already complete"; continue; }
    fi
    echo "  > ${acc}_R${i}"
    curl -fsS -L -C - --retry 5 --retry-delay 3 -o "$dest" "https://${url}" || {
      echo "  !! ${acc}_R${i} download failed" >&2; return 1; }
    if [[ -n "$MD5BIN" ]]; then
      have=$($MD5BIN "$dest" | awk '{print $1}')
      [[ "$have" == "$want" ]] || {
        echo "  !! ${acc}_R${i} MD5 MISMATCH - delete it and rerun" >&2; return 1; }
    fi
  done
  if [[ $DECOMPRESS -eq 1 ]]; then
    for i in 1 2; do
      [[ -s "${OUTDIR}/${acc}_R${i}.fastq" ]] || gunzip -k "${OUTDIR}/${acc}_R${i}.fastq.gz"
    done
  fi
}
export -f fetch_one; export OUTDIR MD5BIN DECOMPRESS

echo "==> downloading with $JOBS parallel jobs"
if command -v parallel >/dev/null; then
  parallel --colsep '\t' -j "$JOBS" fetch_one {1} {2} {3} :::: "$MANIFEST"
else
  running=0
  while IFS=$'\t' read -r acc ftp md5 _; do
    fetch_one "$acc" "$ftp" "$md5" &
    running=$((running+1))
    (( running % JOBS == 0 )) && wait
  done < "$MANIFEST"
  wait
fi

# ---- sample sheet for the pipeline ----
SHEET="samples_new.tsv"
printf 'sample_name\traw_reads_r1\traw_reads_r2\n' > "$SHEET"
ok=0; bad=0
for acc in "${ACCESSIONS[@]}"; do
  if ! r1=$(have_reads "$acc"); then
    echo "  !! $acc missing, left out of $SHEET" >&2; bad=$((bad+1)); continue
  fi
  r2="${r1/_R1./_R2.}"; r2="${r2/_1./_2.}"
  if [[ -s "$r2" ]]; then
    d=$(cd "$(dirname "$r1")" && pwd)
    printf '%s\t%s\t%s\n' "$acc" "$d/$(basename "$r1")" "$d/$(basename "$r2")" >> "$SHEET"
    ok=$((ok+1))
  else
    echo "  !! $acc has R1 but no R2, left out of $SHEET" >&2; bad=$((bad+1))
  fi
done

echo
echo "==> $ok runs ready, $bad incomplete"
echo "    reads : $OUTDIR/"
echo "    sheet : $SHEET"
echo
echo "To extract features, concatenate $SHEET with your existing samples.tsv and"
echo "rerun the pipeline. CNV depth is per-sample against a fixed reference, so"
echo "those columns can simply be appended; the core-SNP alignment has to be"
echo "rebuilt across all samples, because adding genomes redefines the core genome."
