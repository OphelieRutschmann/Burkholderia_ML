# workflow/modules/assembly_illumina.smk
# Workflow for assembly using Illumina reads

rule spades:
    input:
        r1="results/00_QC/fastp/trimmed_reads/{sample}_trimmed_R1.fastq",
        r2="results/00_QC/fastp/trimmed_reads/{sample}_trimmed_R2.fastq"
    output:
        assembly="results/01_Assembly/{sample}/contigs.fasta"
    params:
        outdir="results/01_Assembly/{sample}",
    resources:
        runtime=config["resources"]["general"]["runtime"],
        mem_mb=config["resources"]["general"]["mem_mb"],
        cpus_per_task=config["resources"]["general"]["cpus"]
    container:
        "workflow/containers/spades.sif"
    shell:
        """
        # Perform assembly using spades
        spades.py \
            -o {params.outdir} \
            -1 {input.r1} \
            -2 {input.r2} \
            --threads {resources.cpus_per_task} \

        # Cleanup intermediate files we don't need
        rm -f {params.outdir}/params.txt
        rm -f {params.outdir}/assembly_graph.fastg
        rm -f {params.outdir}/K*-contigs.fasta
        rm -f {params.outdir}/scaffolds.fasta
        rm -f {params.outdir}/scaffolds.paths
        rm -f {params.outdir}/assembly_graph_with_scaffolds.gfa
        rm -rf {params.outdir}/tmp/   
        rm -rf {params.outdir}/K*/
        rm -f {params.outdir}/dataset.info
        rm -rf {params.outdir}/pipeline_state/  
        rm -rf {params.outdir}/misc/  
        rm -f {params.outdir}/assembly_graph_after_simplification.gfa
        rm -f {params.outdir}/aligned_sorted.bam
        rm -f {params.outdir}/aligned_sorted.bam.bai
        rm -f {params.outdir}/aligned.sam
        rm -f {params.outdir}/before_rr.fasta
        rm -f {params.outdir}/run_spades.sh
        rm -f {params.outdir}/run_spades.yaml
        rm -f {params.outdir}/contigs.paths
        rm -f {params.outdir}/coverage.txt
        """

rule filter:
    input:
        assembly="results/01_Assembly/{sample}/contigs.fasta"
    output:
        assembly_filtered="results/01_Assembly/{sample}/contigs_filtered.fasta"
    params:
        outdir="results/01_Assembly/{sample}"
    resources:
        runtime=config["resources"]["general"]["runtime"],
        mem_mb=config["resources"]["general"]["mem_mb"],
        cpus_per_task=config["resources"]["general"]["cpus"]
    container:
        "workflow/containers/seqkit.sif"
    shell:
        """
        # Filter contigs shorter than 500bp using seqkit

        echo "Filtering contigs shorter than 500 bp"
        seqkit seq -m 500 {input.assembly} > {output.assembly_filtered}
        """

rule bwa:
    input:
        r1="results/00_QC/fastp/trimmed_reads/{sample}_trimmed_R1.fastq",
        r2="results/00_QC/fastp/trimmed_reads/{sample}_trimmed_R2.fastq",
        assembly="results/01_Assembly/{sample}/contigs_filtered.fasta"
    output:
        alignement="results/01_Assembly/{sample}/aligned.sam"
    resources:
        runtime=config["resources"]["minimap2"]["runtime"],
        mem_mb=config["resources"]["minimap2"]["mem_mb"],
        cpus_per_task=config["resources"]["minimap2"]["cpus"]
    container:
        "workflow/containers/bwa.sif"
    shell:
        """
        mkdir -p $(dirname {output.alignement})

        # index your assembly
        bwa index {input.assembly}

        bwa mem -t {resources.cpus_per_task} \
            {input.assembly} \
            {input.r1} {input.r2} > {output.alignement}
        """

rule sort_index:
    input:
        alignement="results/01_Assembly/{sample}/aligned.sam"
    output:
        sorted_alignement="results/01_Assembly/{sample}/aligned_sorted.bam",
        bam_index="results/01_Assembly/{sample}/aligned_sorted.bam.bai",
        depth="results/01_Assembly/{sample}/coverage.txt"
    resources:
        runtime=config["resources"]["minimap2"]["runtime"],
        mem_mb=config["resources"]["minimap2"]["mem_mb"],
        cpus_per_task=config["resources"]["minimap2"]["cpus"]
    container:
        "workflow/containers/samtools.sif"
    shell:
        """
        samtools sort -@ {resources.cpus_per_task} {input.alignement} -o {output.sorted_alignement}
        samtools index {output.sorted_alignement}
        """

rule pilon:
    input:
        assembly="results/01_Assembly/{sample}/contigs_filtered.fasta",
        sorted_alignement="results/01_Assembly/{sample}/aligned_sorted.bam",
        bam_index="results/01_Assembly/{sample}/aligned_sorted.bam.bai"
    output:
        polished="results/01_Assembly/{sample}/{sample}_polished.fasta",
        assembly="data/assemblies/{sample}.fasta"
    params:
        outdir="results/01_Assembly/{sample}",
        output_prefix="{sample}_polished"
    resources:
        runtime=config["resources"]["pilon"]["runtime"],
        mem_mb=config["resources"]["pilon"]["mem_mb"],
        cpus_per_task=config["resources"]["pilon"]["cpus"]
    container:
        "workflow/containers/pilon.sif"
    shell:
        """
        mkdir -p {params.outdir}
        mkdir -p $(dirname {output.assembly})

        assembly_path=$(realpath {input.assembly})
        sorted_path=$(realpath {input.sorted_alignement})

        cd {params.outdir}


        pilon --genome "$assembly_path" \
            --frags "$sorted_path" \
            --output {params.output_prefix} \
            --threads {resources.cpus_per_task} \
            --changes \
            --fix all
        cd -

        cp {output.polished} {output.assembly}
        echo "Pilon polishing complete for {wildcards.sample}"
        """