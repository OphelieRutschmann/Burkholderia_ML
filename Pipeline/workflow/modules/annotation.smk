# workflow/modules/annotation.smk
include: "../rules/common.smk"

rule bakta: 
    input:
        assembly="data/assemblies/{sample}.fasta",
        bakta_db="databases/bakta_db/.download_complete"
    output: 
        annotation="results/05_Annotation/{sample}/{sample}.faa"
    params:
        outdir="results/05_Annotation/{sample}",
        db_path=config["databases_bakta"]["bakta_path"]
    resources:
        runtime=config["resources"]["general"]["runtime"], 
        mem_mb=config["resources"]["general"]["mem_mb"],    
        cpus_per_task=config["resources"]["general"]["cpus"] 
    container:
        "workflow/containers/bakta.sif"
    shell:
        """

        # Get absolute path to db
        if [ -n "{params.db_path}" ] && [ -d "{params.db_path}" ]; then
            BAKTA_DB="{params.db_path}"
        else
            BAKTA_DB="databases/bakta_db/db-light"
        fi

        bakta \
            --db $BAKTA_DB \
            --output {params.outdir} \
            --force \
            --prefix {wildcards.sample} \
            {input.assembly}
        """

