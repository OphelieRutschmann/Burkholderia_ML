# workflow/modules/variant_calling_sv.smk
include: "../rules/common.smk"

## This module detects structural variants using Pindel. 
## Before using Pindel, PCR duplicates are marked using Picard

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
        config="results/03_SV/{sample}/pindel_config.txt",
        metrics="results/03_SV/{sample}/insert_metrics.txt"
    params:
        outdir="results/03_SV/{sample}"
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
        # Extract mean insert size: skip comment lines (#) and the column header line
        # then take the first field of the first data line.
        INSERT_SIZE=$(awk 'BEGIN{{FS="\t"}} /^#/{{next}} /^MEDIAN_INSERT_SIZE/{{next}} NF>0{{print int($1); exit}}' {output.metrics})

        echo "insert size: $INSERT_SIZE"

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
        config="results/03_SV/{sample}/pindel_config.txt",
        ref="databases/genomes/refgenome/ref_genome.fasta",
        fai="databases/genomes/refgenome/ref_genome.fasta.fai"
    output:
        flag=touch("results/03_SV/{sample}/.pindel_done")
    params:
        prefix="results/03_SV/{sample}/pindel"
    resources:
        runtime=config["resources"]["general"]["runtime"],
        mem_mb=config["resources"]["general"]["mem_mb"],
        cpus_per_task=config["resources"]["general"]["cpus"]
    container:
        "workflow/containers/pindel.sif"
    shell:
        """
        pindel \
            -f {input.ref} \
            -i {input.config} \
            -o {params.prefix} \
            -T {resources.cpus_per_task}
        """

rule pindel2vcf:
    input:
        ref="databases/genomes/refgenome/ref_genome.fasta",
        flag="results/03_SV/{sample}/.pindel_done",
        pindel="results/03_SV/{sample}/pindel_{svtype}"
    output:
        vcf="results/03_SV/{sample}/cnv_{svtype}.vcf"
    params:
        refname="ref",
        refdate="20240101"
    resources:
        runtime=config["resources"]["general"]["runtime"],
        mem_mb=config["resources"]["general"]["mem_mb"],
        cpus_per_task=config["resources"]["general"]["cpus"]
    container:
        "workflow/containers/pindel.sif"
    shell:
        """
        pindel2vcf \
            -p {input.pindel} \
            -r {input.ref} \
            -R {params.refname} \
            -d {params.refdate} \
            -v {output.vcf}
        """

rule index_cnv_vcf:
    input:
        vcf="results/03_SV/{sample}/cnv_{svtype}.vcf"
    output:
        vcf_gz="results/03_SV/{sample}/cnv_{svtype}.vcf.gz",
        tbi="results/03_SV/{sample}/cnv_{svtype}.vcf.gz.tbi"
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
        vcfs=lambda _: expand(
            "results/03_SV/{sample}/cnv_{svtype}.vcf.gz",
            sample=get_passed_samples(),
            svtype=["TD", "INV", "LI", "D"]
        ),
        tbis=lambda _: expand(
            "results/03_SV/{sample}/cnv_{svtype}.vcf.gz.tbi",
            sample=get_passed_samples(),
            svtype=["TD", "INV", "LI", "D"]
        )
    output:
        merged="results/03_SV/merged_cnv.vcf.gz",
        tbi="results/03_SV/merged_cnv.vcf.gz.tbi"
    resources:
        runtime=config["resources"]["general"]["runtime"],
        mem_mb=config["resources"]["general"]["mem_mb"],
        cpus_per_task=config["resources"]["general"]["cpus"]
    container:
        "workflow/containers/bcftools.sif"
    shell:
        """
        bcftools merge {input.vcfs} -Oz -o {output.merged} --missing-to-ref --force-samples
        tabix -p vcf {output.merged}
        """