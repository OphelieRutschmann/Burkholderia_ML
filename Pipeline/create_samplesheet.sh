#!/bin/bash

reads_dir="/mnt/nfs/bio/Sequencing/Bacteriology/p024_Burkholderia_ML/Genomes"   
output="/home/ruop/Projects/2026_Burkholderia_ML/Pipeline/samples.tsv"

echo -e "sample_name\traw_reads_r1\traw_reads_r2" > "$output"

for r1 in "$reads_dir"/*_R1.fastq; do
    sample=$(basename "$r1" _R1.fastq)
    r2="$reads_dir/${sample}_R2.fastq"

    if [[ -f "$r2" ]]; then
        echo -e "${sample}\t$(realpath "$r1")\t$(realpath "$r2")" >> "$output"
    else
        echo "WARNING: no R2 found for $sample, skipping" >&2
    fi
done

echo "Done. Wrote $(($(wc -l < "$output") - 1)) samples to $output"