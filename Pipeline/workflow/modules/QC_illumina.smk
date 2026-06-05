# workflow/modules/QC_illumina.smk
# Module for QC and trimming of illumina reads

include: "../rules/common.smk"

rule fastQC:
    input: 
        unpack(get_illumina_raw)
    output:
        r1="results/00_QC/fastqc/{sample}_R1_fastqc.html",
        r2="results/00_QC/fastqc/{sample}_R2_fastqc.html",
        r1_zip="results/00_QC/fastqc/{sample}_R1_fastqc.zip",
        r2_zip="results/00_QC/fastqc/{sample}_R2_fastqc.zip"
    params:
        outdir="results/00_QC/fastqc"
    resources:
        runtime=config["resources"]["general"]["runtime"],
        mem_mb=config["resources"]["general"]["mem_mb"],
        cpus_per_task=config["resources"]["general"]["cpus"]
    container:
        "workflow/containers/fastqc.sif"
    shell:
        """
        mkdir -p {params.outdir}
        
        # Create temp directory for FastQC output
        TMPDIR=$(mktemp -d -p {params.outdir})

        fastqc -o $TMPDIR {input.r1} {input.r2}

        mv $TMPDIR/*R1_fastqc.html {output.r1}
        mv $TMPDIR/*R2_fastqc.html {output.r2}
        mv $TMPDIR/*R1_fastqc.zip {output.r1_zip}
        mv $TMPDIR/*R2_fastqc.zip {output.r2_zip}

        rm -rf $TMPDIR
        """

rule fastp:
    input: 
        unpack(get_illumina_raw)
    output:
        html="results/00_QC/fastp/{sample}_fastp.html",
        json="results/00_QC/fastp/{sample}_fastp.json",
        trimmed_r1="results/00_QC/fastp/trimmed_reads/{sample}_trimmed_R1.fastq",
        trimmed_r2="results/00_QC/fastp/trimmed_reads/{sample}_trimmed_R2.fastq"
    params:
        outdir="results/00_QC/fastp/trimmed_reads",
    container:
        "workflow/containers/fastp.sif"
    resources:
        runtime=config["resources"]["general"]["runtime"],
        mem_mb=config["resources"]["general"]["mem_mb"],
        cpus_per_task=config["resources"]["general"]["cpus"]
    shell:
        """
        mkdir -p {params.outdir}
        
        fastp \
            --in1 {input.r1} \
            --in2 {input.r2} \
            --out1 {output.trimmed_r1} \
            --out2 {output.trimmed_r2} \
            --length_required 100 \
            --detect_adapter_for_pe \
            --qualified_quality_phred 30 \
            --html {output.html} \
            --json {output.json}
        """   

rule multiqc:
    input:
        expand("results/00_QC/fastp/{sample}_fastp.json", sample=config["samples"]["name"])
    output:
        report="results/00_QC/multiqc/multiqc_report.html",
        data="results/00_QC/multiqc/multiqc_data/multiqc_fastp.txt"
    resources:
        runtime=config["resources"]["general"]["runtime"],
        mem_mb=config["resources"]["general"]["mem_mb"],
        cpus_per_task=config["resources"]["general"]["cpus"]
    container:
        "workflow/containers/multiqc.sif"
    shell:
        """
        multiqc results/00_QC/fastp/ -o results/00_QC/multiqc/ --force
        """

rule bwa_index:
    input:
        ref="databases/genomes/refgenome/ref_genome.fasta"
    output:
        amb="databases/genomes/refgenome/ref_genome.fasta.amb",
        ann="databases/genomes/refgenome/ref_genome.fasta.ann",
        bwt="databases/genomes/refgenome/ref_genome.fasta.bwt",
        pac="databases/genomes/refgenome/ref_genome.fasta.pac",
        sa="databases/genomes/refgenome/ref_genome.fasta.sa"
    resources:
        runtime=config["resources"]["general"]["runtime"],
        mem_mb=config["resources"]["general"]["mem_mb"],
        cpus_per_task=config["resources"]["general"]["cpus"]
    container:
        "workflow/containers/bwa.sif"
    shell:
        """
        bwa index {input.ref}
        """

# Align to reference genome to calculate coverage
rule bwa:
    input:
        r1="results/00_QC/fastp/trimmed_reads/{sample}_trimmed_R1.fastq",
        r2="results/00_QC/fastp/trimmed_reads/{sample}_trimmed_R2.fastq",
        ref="databases/genomes/refgenome/ref_genome.fasta",
        index="databases/genomes/refgenome/ref_genome.fasta.ann",
    output:
        alignement="results/01_alignement/{sample}.aligned.sam"
    resources:
        runtime=config["resources"]["general"]["runtime"],
        mem_mb=config["resources"]["general"]["mem_mb"],
        cpus_per_task=config["resources"]["general"]["cpus"]
    container:
        "workflow/containers/bwa.sif"
    shell:
        """
        mkdir -p $(dirname {output.alignement})

        bwa mem -t {resources.cpus_per_task} \
            {input.ref} \
            {input.r1} {input.r2} > {output.alignement}
        """

rule coverage:
    input:
        alignement="results/01_alignement/{sample}.aligned.sam"
    output:
        bam="results/01_alignement/{sample}.aligned.sorted.bam",      
        bai="results/01_alignement/{sample}.aligned.sorted.bam.bai",
        coverage="results/00_QC/coverage/{sample}/coverage.txt",
        mapping_stats="results/00_QC/coverage/{sample}/mapping_stats.txt"
    params:
        outdir="results/00_QC/coverage/{sample}",
        alignement_dir="results/01_alignement"
    resources:
        runtime=config["resources"]["general"]["runtime"],
        mem_mb=config["resources"]["general"]["mem_mb"],
        cpus_per_task=config["resources"]["general"]["cpus"]
    container:
        "workflow/containers/samtools.sif"
    shell:
        """
        mkdir -p {params.outdir}

        samtools view -bS {input.alignement} | samtools sort -o {output.bam}
        samtools index {output.bam}

        samtools coverage {output.bam} -o {output.coverage}
        samtools flagstat {output.bam} > {output.mapping_stats}
        """

checkpoint filter:
    input: 
        multiqc="results/00_QC/multiqc/multiqc_data/multiqc_fastp.txt",
        coverage=expand("results/00_QC/coverage/{sample}/coverage.txt", sample=config["samples"]["name"])
    output:
        "results/00_QC/passed_samples.csv"
    resources:
        runtime=config["resources"]["general"]["runtime"],
        mem_mb=config["resources"]["general"]["mem_mb"],
        cpus_per_task=config["resources"]["general"]["cpus"]
    script:
        "../scripts/filter_reads.py"