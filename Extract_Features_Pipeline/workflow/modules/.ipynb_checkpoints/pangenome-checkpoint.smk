# workflow/modules/pangenome.smk
# Pangenome analysis: gene presence/absence across all strains using Roary

include: "../rules/common.smk"

rule annotate_strain:
    input:
        genome=config.get("genomes_dir", "databases/genomes") + "/" +
               config.get("genome_pattern", "{sample}/{sample}.fasta"),
        flag="databases/bakta_db/.download_complete"
    output:
        gbff="databases/annotations/{sample}.gbff"
    params:
        outdir="databases/annotations/{sample}",
        db_path=config["databases_bakta"]["bakta_path"]
    threads:
        config["resources"]["general"]["cpus"]
    resources:
        runtime=config["resources"]["general"]["runtime"],
        mem_mb=config["resources"]["general"]["mem_mb"],
        cpus_per_task=config["resources"]["general"]["cpus"]
    container:
        "workflow/containers/bakta.sif"
    shell:
        """
        mkdir -p {params.outdir}

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
            --keep-contig-headers \
            {input.genome}

        # Move GBFF to expected location
        mv {params.outdir}/{wildcards.sample}.gbff {output.gbff}
        """

rule strain_gbff_to_gff3:
    input:
        gbff="databases/annotations/{sample}.gbff"
    output:
        gff3=config.get("annotation_dir", "databases/annotations") + "/" +
             config.get("annotation_pattern", "{sample}.gff")
    container:
        "workflow/containers/biopython.sif"
    resources:
        runtime=config["resources"]["general"]["runtime"],
        mem_mb=config["resources"]["general"]["mem_mb"],
        cpus_per_task=config["resources"]["general"]["cpus"]
    shell:
        """
        python workflow/scripts/gbff_gff3.py \
            --input {input.gbff} \
            --output {output.gff3}
        """

rule run_roary:
    input:
        gff_files=expand(
            config.get("annotation_dir", "databases/annotations") + "/" +
            config.get("annotation_pattern", "{sample}.gff"),
            sample=get_passed_samples()
        )
    output:
        gene_presence_absence="results/05_Pangenome/gene_presence_absence.csv",
        summary_statistics="results/05_Pangenome/summary_statistics.txt"
    params:
        outdir="results/05_Pangenome",
        identity=config.get("roary", {}).get("identity", 95),
        cd=config.get("roary", {}).get("cd", 95)
    threads:
        config["resources"]["general"]["cpus"]
    resources:
        runtime=config["resources"]["general"]["runtime"],
        mem_mb=config["resources"]["general"]["mem_mb"],
        cpus_per_task=config["resources"]["general"]["cpus"]
    container:
        "workflow/containers/roary.sif"
    shell:
        """
        mkdir -p {params.outdir}

        roary -f {params.outdir} \
              -e -n -v \
              -i {params.identity} \
              -cd {params.cd} \
              -t {threads} \
              {input.gff_files}

        echo "Pangenome analysis complete"
        echo "Genes found: $(grep -v '^Gene' {output.gene_presence_absence} | wc -l)"
        """

rule pangenome_to_feature_matrix:
    input:
        gene_matrix="results/05_Pangenome/gene_presence_absence.csv"
    output:
        feature_matrix="results/genes_presence_absence.tsv"
    container:
        "workflow/containers/python.sif"
    shell:
        """
        python workflow/scripts/convert_pangenome_to_matrix.py \
            --input {input.gene_matrix} \
            --output {output.feature_matrix}
        """