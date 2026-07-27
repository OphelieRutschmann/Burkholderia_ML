# workflow/modules/databases.smk
include: "../rules/common.smk"

# Download the bakta database 
rule download_bakta:
    output:
        flag="databases/bakta_db/.download_complete"
    params:
        db_type="light",
        db_path="databases/bakta_db",
        custom_path=config["databases_bakta"]["bakta_path"]
    container:
        "workflow/containers/bakta.sif"
    resources:
        runtime=config["resources"]["general"]["runtime"],
        mem_mb=config["resources"]["general"]["mem_mb"],
        cpus_per_task=config["resources"]["general"]["cpus"]
    shell:
        """
        # Use provided bakta database if it exists
        if [ -d "{params.db_path}" ] && [ "$(ls -A {params.db_path})" ]; then
            echo "Using provided Bakta database"
            touch {output.flag}
            exit 0
        fi

        # Else download the Bakta database
        bakta_db download --output {params.db_path} --type {params.db_type}
        
        # Create completion flag
        touch {output.flag}
        """