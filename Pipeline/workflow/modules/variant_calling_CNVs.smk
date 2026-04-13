# workflow/modules/variant_calling_snps.smk
include: "../rules/common.smk"

## This module detects CNV (copy number variations) using Pindel. Before using Pindel, PCR duplicates are marked using Picard ##
rule sort_bam:
  input:
    bam="results/mapping/{sample}.bam"
  output:
    bam="results/mapping/{sample}.sorted.bam",
    bai="results/mapping/{sample}.sorted.bam.bai"
  container:
    "workflow/containers/samtools.sif"
  shell:
    """
    samtools sort -o {output.bam} {input.bam}
    samtools index {output.bam}
    """

rule mark_duplicates:
  input:
    bam="results/mapping/{sample}.sorted.bam"
  output:
    bam="results/mapping/{sample}.dedup.bam",
    bai="results/mapping/{sample}.dedup.bam.bai",
    metrics="results/mapping/{sample}.dedup.metrics.txt"
  container:
    "workflow/containers/picard.sif"
  shell:
    """
    picard MarkDuplicates \
      I={input.bam} \
      O={output.bam} \
      M={output.metrics} \
      REMOVE_DUPLICATES=false \
      CREATE_INDEX=true \
      VALIDATION_STRINGENCY=SILENT
      """

rule pindel_config:
  input:
    bam="results/mapping/{sample}.dedup.bam"
  output:
    config="results/CNV/{sample}/pindel_config.txt",
    metrics="results/CNV/{sample}/insert_metrics.txt"
  params:
    outdir="results/CNV/{sample}"
  container:
    "workflow/containers/picard.sif"
  shell:
    """
    mkdir -p {params.outdir}
        
    # Compute insert size metrics
    picard CollectInsertSizeMetrics \
      I={input.bam} \
      O={output.metrics} \
      H={params.outdir}/insert_histogram.pdf

    # Extract mean insert size from the insert_metrics.txt file (first data line after header)
    INSERT_SIZE=$(awk 'BEGIN{{FS="\\t"}} !/^#/ && NR==8 {{print int($1)}}' {output.metrics})

    echo "insert size: $INSERT_SIZE"

    # Create pindel config file
    echo "{input.bam} $INSERT_SIZE {wildcards.sample}" > {output.config}
    """

rule pindel:
    input:
      config="results/CNV/{sample}/pindel_config.txt",
      ref="databases/genomes/refgenome/ref_genome.fasta"
    output:
      prefix="results/CNV/{sample}/pindel"
    container:
      "workflow/containers/pindel.sif"
    threads: 4
    shell:
      """
      pindel \
        -f {input.ref} \
        -i {input.config} \
        -o {output.prefix} \
        -T {threads}
      """

rule pindel_to_vcf_cnv:
  input:
    prefix="results/CNV/{sample}/pindel",
    ref="databases/genomes/refgenome/ref_genome.fasta"
  output:
    vcf="results/CNV/{sample}/cnv.vcf"
  container:
    "workflow/containers/pindel.sif"
  shell:
    """
    pindel2vcf \
      -p {input.prefix}_TD \
      -r {input.ref} \
      -R ref \
      -d 20240101 \
      -v {output.vcf}
    """

rule index_cnv_vcf:
  input:
    vcf="results/CNV/{sample}/cnv.vcf"
  output:
    vcf_gz="results/CNV/{sample}/cnv.vcf.gz",
    tbi="results/CNV/{sample}/cnv.vcf.gz.tbi"
  container:
    "workflow/containers/bcftools.sif"
  shell:
    """
    bcftools view -Oz -o {output.vcf_gz} {input.vcf}
    tabix -p vcf {output.vcf_gz}
    """

rule merge_cnv_vcfs:
  input:
    expand("results/CNV/{sample}/cnv.vcf.gz", sample=SAMPLES)
  output:
    merged="results/CNV/merged_cnv.vcf.gz"
  container:
    "workflow/containers/bcftools.sif"
  shell:
    """
    bcftools merge {input} -Oz -o {output.merged}
    tabix -p vcf {output.merged}
    """