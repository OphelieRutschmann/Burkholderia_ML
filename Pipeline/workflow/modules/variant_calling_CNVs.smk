# workflow/modules/variant_calling_cnv.smk
include: "../rules/common.smk"

## This module detects CNV (copy number variations) using Pindel. Before using Pindel, PCR duplicates are marked using Picard

# Use alignement to reference generated in the QC_Illumina.smk module (for coverage calculation)
rule mark_duplicates:
  input:
    bam="results/alignement/{sample}.aligned.sorted.bam"
  output:
    bam="results/alignement/{sample}.dedup.bam",
    bai="results/alignement/{sample}.dedup.bam.bai",
    metrics="results/alignement/{sample}.dedup.metrics.txt"
  resources:
    runtime=config["resources"]["general"]["runtime"],
    mem_mb=config["resources"]["general"]["mem_mb"],
    cpus_per_task=config["resources"]["general"]["cpus"]
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
    bam="results/alignement/{sample}.dedup.bam"
  output:
    config="results/CNV/{sample}/pindel_config.txt",
    metrics="results/CNV/{sample}/insert_metrics.txt"
  params:
    outdir="results/CNV/{sample}"
  resources:
    runtime=config["resources"]["general"]["runtime"],
    mem_mb=config["resources"]["general"]["mem_mb"],
    cpus_per_task=config["resources"]["general"]["cpus"]
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
    resources:
      runtime=config["resources"]["general"]["runtime"],
      mem_mb=config["resources"]["general"]["mem_mb"],
      cpus_per_task=config["resources"]["general"]["cpus"]
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
  resources:
    runtime=config["resources"]["general"]["runtime"],
    mem_mb=config["resources"]["general"]["mem_mb"],
    cpus_per_task=config["resources"]["general"]["cpus"]
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
  resources:
    runtime=config["resources"]["general"]["runtime"],
    mem_mb=config["resources"]["general"]["mem_mb"],
    cpus_per_task=config["resources"]["general"]["cpus"]
  container:
    "workflow/containers/bcftools.sif"
  shell:
    """
    bcftools view -Oz -o {output.vcf_gz} {input.vcf}
    tabix -p vcf {output.vcf_gz}
    """

rule merge_cnv_vcfs:
  input:
    lambda _: expand(
      "results/CNV/{sample}/cnv.vcf.gz",
      sample=get_passed_samples()
    )
  output:
    merged="results/CNV/merged_cnv.vcf.gz"
  resources:
    runtime=config["resources"]["general"]["runtime"],
    mem_mb=config["resources"]["general"]["mem_mb"],
    cpus_per_task=config["resources"]["general"]["cpus"]
  container:
    "workflow/containers/bcftools.sif"
  shell:
    """
    bcftools merge {input} -Oz -o {output.merged}
    tabix -p vcf {output.merged}
    """