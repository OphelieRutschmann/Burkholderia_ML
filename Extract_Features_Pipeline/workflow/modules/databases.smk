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
# Copy the reference genome given in config/databases.yaml (ref_genome) into
# databases/genomes/refgenome/. These rules only run when the file is missing there.
rule ref_genome_fasta:
    output:
        fasta="databases/genomes/refgenome/ref_genome.fasta"
    params:
        src=config.get("ref_genome", {}).get("fasta", "")
    resources:
        runtime=res("light", "runtime"),
        mem_mb=res("light", "mem_mb"),
        cpus_per_task=res("light", "cpus")
    shell:
        """
        if [ -z "{params.src}" ] || [ ! -f "{params.src}" ]; then
            echo "Reference genome not found. Set ref_genome: fasta: in config/databases.yaml" >&2
            echo "to the reference FASTA, or put it at {output.fasta}" >&2
            exit 1
        fi
        cp "{params.src}" {output.fasta}
        """

rule ref_genome_gbff:
    output:
        gbff="databases/genomes/refgenome/ref_genome.gbff"
    params:
        src=config.get("ref_genome", {}).get("gbff", "")
    resources:
        runtime=res("light", "runtime"),
        mem_mb=res("light", "mem_mb"),
        cpus_per_task=res("light", "cpus")
    shell:
        """
        if [ -z "{params.src}" ] || [ ! -f "{params.src}" ]; then
            echo "Reference annotation not found. Set ref_genome: gbff: in config/databases.yaml" >&2
            echo "to the reference GenBank file, or put it at {output.gbff}" >&2
            exit 1
        fi
        cp "{params.src}" {output.gbff}
        """
