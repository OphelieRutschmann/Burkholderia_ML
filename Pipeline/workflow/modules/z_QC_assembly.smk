# workflow/modules/QC_assembly.smk
# Module for assembly quality assessment

# Include shared functions
include: "../rules/common.smk"

rule busco:
    input:
        assembly="data/assemblies/{sample}.fasta",
        flag="databases/busco/.download_complete"
    output:
        txt="results/02_Assembly_QC/02a_Busco/{sample}/short_summary.{sample}.txt",
        json="results/02_Assembly_QC/02a_Busco/{sample}/short_summary.{sample}.json"
    resources:
        runtime=config["resources"]["busco"]["runtime"],
        mem_mb=config["resources"]["busco"]["mem_mb"],
        cpus_per_task=config["resources"]["busco"]["cpus"]
    params:
        outdir="results/02_Assembly_QC/02a_Busco/{sample}",
        busco_summaries="results/02_Assembly_QC/02a_Busco/Busco_summaries",
        lineage=config["quality"]["busco"]["lineage"],
        db_path=config["databases_busco"]["busco_path"]
    container:
        "workflow/containers/busco.sif"
    shell:
        """
        # Create output directories
        mkdir -p {params.outdir}
        mkdir -p {params.busco_summaries}

        # Get absolute path to input file
        INPUT_ABS=$(realpath {input.assembly})

        # Get absolute path to db
        if [ -n "{params.db_path}" ] && [ -d "{params.db_path}" ]; then
            BUSCO_LINEAGE_PATH="{params.db_path}"
        else
            BUSCO_LINEAGE_PATH=$(realpath databases/busco/{params.lineage})
        fi
        
        # Run BUSCO
        cd {params.outdir}

        busco -i $INPUT_ABS \
            -m genome \
            --lineage_dataset ${{BUSCO_LINEAGE_PATH}} \
            -o busco_run \
            -f \
            -c {resources.cpus_per_task} \
            --offline
        
        # Rename the summary file to match expected output
        mv busco_run/run_{params.lineage}/short_summary.txt short_summary.{wildcards.sample}.txt
        mv busco_run/run_{params.lineage}/short_summary.json short_summary.{wildcards.sample}.json

        # Copy summary file to summaries directory
        cp short_summary.{wildcards.sample}.txt ../Busco_summaries/
        cp short_summary.{wildcards.sample}.json ../Busco_summaries/
        """

rule busco_summaries:
    input:
        # ensures busco_summaries waits for all busco jobs to complete
        expand("results/02_Assembly_QC/02a_Busco/{sample}/short_summary.{sample}.json", 
               sample=config["samples"]["name"]) 
    output:
        report(
            "results/02_Assembly_QC/02a_Busco/Busco_summaries/busco_figure.png",
            caption="report/busco.rst",
            category="QC assembly"
        )
    params: 
        busco_dir="results/02_Assembly_QC/02a_Busco/Busco_summaries"
    resources:
        runtime=config["resources"]["general"]["runtime"],
        mem_mb=config["resources"]["general"]["mem_mb"],
        cpus_per_task=config["resources"]["general"]["cpus"]
    container:
        "workflow/containers/busco.sif"
    shell:
        """
        busco --plot {params.busco_dir}
        """

rule quast:
    input:
        assembly="data/assemblies/{sample}.fasta"
    output:
        "results/02_Assembly_QC/02b_Quast/{sample}/report.txt"
    resources:
        runtime=config["resources"]["general"]["runtime"],
        mem_mb=config["resources"]["general"]["mem_mb"],
        cpus_per_task=config["resources"]["general"]["cpus"]
    params:
        outdir="results/02_Assembly_QC/02b_Quast"
    container:
        "workflow/containers/quast.sif"
    shell:
        """      
        quast.py \
            -o {params.outdir}/{wildcards.sample} \
            --gene-finding \
            -t {resources.cpus_per_task} \
            {input.assembly}

        rm -f results/02_Assembly_QC/02b_Quast/{wildcards.sample}/transposed_report*
        rm -f results/02_Assembly_QC/02b_Quast/{wildcards.sample}/icarus.html
        rm -f results/02_Assembly_QC/02b_Quast/{wildcards.sample}/report.tex
        rm -rf results/02_Assembly_QC/02b_Quast/{wildcards.sample}/icarus_viewers/
        rm -rf results/02_Assembly_QC/02b_Quast/{wildcards.sample}/predicted_genes/
        rm -rf results/02_Assembly_QC/02b_Quast/{wildcards.sample}/basic_stats/
        """

rule CheckM2:
    input: 
        assembly="data/assemblies/{sample}.fasta",
        checkm_flag="databases/CheckM2_database/.download_complete"
    output:
        "results/02_Assembly_QC/02c_CheckM/{sample}/quality_report.tsv"
    params:
        db_path=config["databases_checkM2"]["checkM2_path"],
        outdir="results/02_Assembly_QC/02c_CheckM/{sample}"
    resources:
        runtime=config["resources"]["checkM"]["runtime"],
        mem_mb=config["resources"]["checkM"]["mem_mb"],
        cpus_per_task=config["resources"]["checkM"]["cpus"]
    container:
        "workflow/containers/checkm2.sif"
    shell:
        """
        # Get absolute path to db
        if [ -n "{params.db_path}" ] && [ -f "{params.db_path}" ]; then
            CHECKM2_PATH="{params.db_path}"
        else
            CHECKM2_PATH="databases/CheckM2_database/checkm.dmnd"
        fi

        echo "$CHECKM2_PATH"
        
        checkm2 predict \
            -i {input.assembly} \
            -o {params.outdir} \
            --threads {resources.cpus_per_task} \
            --database_path $CHECKM2_PATH \
            --force
        """


checkpoint filter_assemblies:
    input:
        checkm =expand("results/02_Assembly_QC/02c_CheckM/{sample}/quality_report.tsv", sample=config["samples"]["name"]),
        quast=expand("results/02_Assembly_QC/02b_Quast/{sample}/report.txt", sample=config["samples"]["name"]),
        busco=expand("results/02_Assembly_QC/02a_Busco/{sample}/short_summary.{sample}.json", sample=config["samples"]["name"])
    output:
        passed_samples="results/02_Assembly_QC/passed_samples.txt",
        report="results/02_Assembly_QC/assembly_qc_report.tsv"
    params:
        samples=config["samples"]["name"],
    resources:
        runtime=config["resources"]["general"]["runtime"],
        mem_mb=config["resources"]["general"]["mem_mb"],
        cpus_per_task=config["resources"]["general"]["cpus"]
    script:
        "../scripts/filter_assemblies.py"