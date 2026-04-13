# workflow/modules/QC_illumina.smk
# Module for QC and trimming of illumina reads

rule fastQC:
    input: 
        r1="results/00_QC/fastp/trimmed_reads/{sample}_trimmed_R1.fastq",
        r2="results/00_QC/fastp/trimmed_reads/{sample}_trimmed_R2.fastq"
    output:
        r1="results/00_QC/fastqc/{sample}_R1_fastqc.html",
        r2="results/00_QC/fastqc/{sample}_R2_fastqc.html"
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

        rm -rf $TMPDIR
        """

rule fastp:
    input: 
        r1="results/00_QC/fastp/trimmed_reads/{sample}_trimmed_R1.fastq",
        r2="results/00_QC/fastp/trimmed_reads/{sample}_trimmed_R2.fastq"
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

rule bwa_index:
    input:
        ref="databases/genomes/refgenome/ref_genome.fasta"
    output:
        amb="databases/genomes/refgenome/ref_genome.fasta.amb",
        ann="databases/genomes/refgenome/ref_genome.fasta.ann",
        bwt="databases/genomes/refgenome/ref_genome.fasta.bwt",
        pac="databases/genomes/refgenome/ref_genome.fasta.pac",
        sa="databases/genomes/refgenome/ref_genome.fasta.sa"
    container:
        "workflow/containers/bwa.sif"
    shell:
        """
        bwa index {input.ref}
        """

rule bwa:
    input:
        r1="results/00_QC/fastp/trimmed_reads/{sample}_trimmed_R1.fastq",
        r2="results/00_QC/fastp/trimmed_reads/{sample}_trimmed_R2.fastq",
        ref="databases/genomes/refgenome/ref_genome.fasta",
        index="databases/genomes/refgenome/ref_genome.fasta.ann",
    output:
        alignement="results/00_QC/coverage/{sample}/aligned.sam"
    resources:
        runtime=config["resources"]["minimap2"]["runtime"],
        mem_mb=config["resources"]["minimap2"]["mem_mb"],
        cpus_per_task=config["resources"]["minimap2"]["cpus"]
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
        alignement="results/00_QC/coverage/{sample}/aligned.sam"
    output:
        coverage="results/00_QC/coverage/{sample}/coverage.txt"
        mapping_stats="results/00_QC/coverage/{sample}/mapping_stats.txt"
    params:
        outdir="results/00_QC/coverage/{sample}"
    resources:
        runtime=config["resources"]["minimap2"]["runtime"],
        mem_mb=config["resources"]["minimap2"]["mem_mb"],
        cpus_per_task=config["resources"]["minimap2"]["cpus"]
    container:
        "workflow/containers/samtools.sif"
    shell:
        """
        mkdir -p {params.outdir}
        
        # Sort and Index the alignement
        samtools view -bS {input.alignement} | samtools sort -o {params.outdir}/output.sorted.bam 
        samtools index {params.outdir}/output.sorted.bam

        # Calculate coverage
        samtools coverage {params.outdir}/output.sorted.bam -o {output.coverage}
        samtools flagstat {params.outdir}/output.sorted.bam > {output.mapping_stats}
        """

checkpoint filter:
    input:
    output:
    params:
        samples=config["samples"]["name"],
    resources:
        runtime=config["resources"]["general"]["runtime"],
        mem_mb=config["resources"]["general"]["mem_mb"],
        cpus_per_task=config["resources"]["general"]["cpus"]
    script:
        "../scripts/filter_reads.py"