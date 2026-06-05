# workflow/modules/variant_calling_cnv.smk
include: "../rules/common.smk"

## This module detects 02_CNV (copy number variations) using Pindel. Before using Pindel, PCR duplicates are marked using Picard

# Use alignement to reference generated in the QC_Illumina.smk module (for coverage calculation)
rule mark_duplicates:
    input:
        bam="results/01_alignement/{sample}.aligned.sorted.bam"
    output:
        bam="results/01_alignement/{sample}.dedup.bam",
        bai="results/01_alignement/{sample}.dedup.bai",
        metrics="results/01_alignement/{sample}.dedup.metrics.txt"
    resources:
        runtime=config["resources"]["general"]["runtime"],
        mem_mb=config["resources"]["general"]["mem_mb"],
        cpus_per_task=config["resources"]["general"]["cpus"]
    container:
        "workflow/containers/picard.sif"
    shell:
            """
            java -jar /usr/picard/picard.jar MarkDuplicates \
                I={input.bam} \
                O={output.bam} \
                M={output.metrics} \
                REMOVE_DUPLICATES=false \
                CREATE_INDEX=true \
                VALIDATION_STRINGENCY=SILENT
            """

rule pindel_config:
    input:
        bam="results/01_alignement/{sample}.dedup.bam"
    output:
        config="results/03_CNV/{sample}/pindel_config.txt",
        metrics="results/03_CNV/{sample}/insert_metrics.txt"
    params:
        outdir="results/03_CNV/{sample}"
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
        java -jar /usr/picard/picard.jar CollectInsertSizeMetrics \
            I={input.bam} \
            O={output.metrics} \
            H={params.outdir}/insert_histogram.pdf

        # Extract mean insert size from the insert_metrics.txt file (first data line after header)
        INSERT_SIZE=$(awk 'BEGIN{{FS="\\t"}} !/^#/ && NR==8 {{print int($1)}}' {output.metrics})

        echo "insert size: $INSERT_SIZE"

        # Create pindel config file
        echo "{input.bam} $INSERT_SIZE {wildcards.sample}" > {output.config}
        """
rule fasta_index:
    input:
        ref="databases/genomes/refgenome/ref_genome.fasta"
    output:
        fai="databases/genomes/refgenome/ref_genome.fasta.fai"
    resources:
        runtime=config["resources"]["general"]["runtime"],
        mem_mb=config["resources"]["general"]["mem_mb"],
        cpus_per_task=config["resources"]["general"]["cpus"]
    container:
        "workflow/containers/samtools.sif"
    shell:
        """
        samtools faidx {input.ref}
        """

rule pindel:
    input:
        config="results/03_CNV/{sample}/pindel_config.txt",
        ref="databases/genomes/refgenome/ref_genome.fasta",
        fai="databases/genomes/refgenome/ref_genome.fasta.fai"
    output:
        flag=touch("results/03_CNV/{sample}/.pindel_done")
    params:
        prefix="results/03_CNV/{sample}/pindel"
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
            -o {params.prefix} \
            -T {threads}
        """

rule pindel_to_vcf_cnv:
    input:
        flag="results/03_CNV/{sample}/.pindel_done",
        ref="databases/genomes/refgenome/ref_genome.fasta"
    output:
        vcf="results/03_CNV/{sample}/cnv.vcf"
    params:
        prefix="results/03_CNV/{sample}/pindel"
    resources:
        runtime=config["resources"]["general"]["runtime"],
        mem_mb=config["resources"]["general"]["mem_mb"],
        cpus_per_task=config["resources"]["general"]["cpus"]
    container:
        "workflow/containers/pindel.sif"
    shell:
        """
        pindel2vcf \
            -p {params.prefix}_TD \
            -r {input.ref} \
            -R ref \
            -d 20240101 \
            -v {output.vcf}
        """

rule index_cnv_vcf:
    input:
        vcf="results/03_CNV/{sample}/cnv.vcf"
    output:
        vcf_gz="results/03_CNV/{sample}/cnv.vcf.gz",
        tbi="results/03_CNV/{sample}/cnv.vcf.gz.tbi"
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
            "results/03_CNV/{sample}/cnv.vcf.gz",
            sample=get_passed_samples()
        )
    output:
        merged="results/03_CNV/merged_cnv.vcf.gz"
    resources:
        runtime=config["resources"]["general"]["runtime"],
        mem_mb=config["resources"]["general"]["mem_mb"],
        cpus_per_task=config["resources"]["general"]["cpus"]
    container:
        "workflow/containers/bcftools.sif"
    shell:
        """
        bcftools merge {input} -Oz -o {output.merged} --missing-to-ref --force-single
        tabix -p vcf {output.merged}
        """