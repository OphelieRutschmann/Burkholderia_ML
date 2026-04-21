# workflow/modules/databases.smk
include: "../rules/common.smk"

# Download reference genome that will be used for variant calling
rule download_reference:
    output:
        ref="databases/genomes/refgenome/ref_genome.fasta",
        ref_gff="databases/genomes/refgenome/ref_genome.gff"
    params:
        ref_gen_dir="databases/genomes/refgenome",
        accession=config["ref_mapping"]["accession"],
        path=config["ref_mapping"]["path"]
    resources:
        runtime=config["resources"]["general"]["runtime"],
        mem_mb=config["resources"]["general"]["mem_mb"],
        cpus_per_task=config["resources"]["general"]["cpus"]
    container:
        "workflow/containers/ncbi.sif"
    shell:
        """
        if [ -f "{output.ref}" ]; then
           echo "Reference genome already exists, skipping download."
           exit 0
        fi

        # If custom database path provided, copy to output dir
        if [ -n "{params.path}" ] && [ -d "{params.path}" ]; then
            echo "Using custom reference from {params.path} for variant calling"
            
            # Copy files to output dir
            cp {params.path} {output.ref}/

            echo "Reference for core SNP setup complete."
            exit 0
        fi

        # Download reference genome
        echo "Downloading reference genome"
        ncbi-genome-download bacteria \
            -F fasta \
            --assembly-accessions "{params.accession}" \
            --flat-output \
            -o {params.ref_gen_dir}
        
        ncbi-genome-download bacteria \
            -F gff \
            --assembly-accessions "{params.accession}" \
            --flat-output \
            -o {params.ref_gen_dir}
                
        # Uncompress and rename the files
        find {params.ref_gen_dir} -name "*.fna.gz" -exec gunzip -c {{}} \\; > {output.ref}
        find {params.ref_gen_dir} -name "*.gff" -exec cat {{}} \; > {output.ref_gff}
        """

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