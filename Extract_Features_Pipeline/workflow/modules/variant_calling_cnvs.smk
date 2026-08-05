# workflow/modules/variant_calling_cnv_depth.smk

include: "../rules/common.smk"

rule gbff_to_gff3:
    input:
        gbff="databases/genomes/refgenome/ref_genome.gbff"
    output:
        gff3="databases/genomes/refgenome/ref_genome.gff3"
    container:
        "workflow/containers/python.sif"
    resources:
        runtime=config["resources"]["general"]["runtime"],
        mem_mb=config["resources"]["general"]["mem_mb"],
        cpus_per_task=config["resources"]["general"]["cpus"]
    shell:
        """
        python workflow/scripts/gbff_to_gff3.py \
            --input {input.gbff} \
            --output {output.gff3}
        """

rule prepare_regions:
    input:
        gff="databases/genomes/refgenome/ref_genome.gff3",
        fai="databases/genomes/refgenome/ref_genome.fasta.fai"
    output:
        bed="databases/genomes/refgenome/regions.bed"
    params:
        promoter_window=config["cnv_depth"].get("promoter_window", 200)
    resources:
        runtime=config["resources"]["general"]["runtime"],
        mem_mb=config["resources"]["general"]["mem_mb"],
        cpus_per_task=config["resources"]["general"]["cpus"]
    container:
        "workflow/containers/bedtools.sif"
    shell:
        """
        # Extract all genomic features from the GFF3 file.
        # Strips comment lines, then parses the 9th attribute column to recover
        # the feature ID and Name, building a label of the form ID__Name__type.
        grep -v "^#" {input.gff} | \
        awk 'BEGIN{{OFS="\t"}} {{

            id="unknown"
            name="unknown"

            n=split($9,a,";")

            for(i=1;i<=n;i++) {{
                split(a[i],kv,"=")

                if(length(kv[2]) > 0) {{
                    if(id=="unknown") id=kv[2]
                    name=kv[2]
                }}
            }}

            feature = id"__"name"__"$3
            print $1,$4-1,$5,feature,".",$7
        }}' \
        > {output.bed}.tmp_features

        # Extract promoter regions for CDS and gene features only.
        # bedtools flank extends each feature upstream to approximate the promoter window.
        grep -v "^#" {input.gff} | \
        awk '$3=="CDS" || $3=="gene"' | \
        awk 'BEGIN{{OFS="\t"}} {{

            id="unknown"
            name="unknown"

            n=split($9,a,";")

            for(i=1;i<=n;i++) {{
                split(a[i],kv,"=")

                if(length(kv[2]) > 0) {{
                    if(id=="unknown") id=kv[2]
                    name=kv[2]
                }}
            }}

            feature = id"__"name"__"$3
            print $1,$4-1,$5,feature,".",$7
        }}' | \
        bedtools flank \
            -i stdin \
            -g {input.fai} \
            -l {params.promoter_window} \
            -r 0 \
            -s | \
        awk 'BEGIN{{OFS="\t"}} {{
            print $1,$2,$3,$4"__promoter",$5,$6
        }}' \
        > {output.bed}.tmp_promoters

        # Merge features and promoters, sort, and deduplicate.
        cat {output.bed}.tmp_features {output.bed}.tmp_promoters | \
        bedtools sort -i stdin | \
        uniq > {output.bed}

        rm {output.bed}.tmp_features {output.bed}.tmp_promoters
        """

rule mosdepth_regions:
    input:
        bam="results/01_alignement/{sample}.dedup.bam",
        bai="results/01_alignement/{sample}.dedup.bai",
        bed="databases/genomes/refgenome/regions.bed"
    output:
        regions="results/04_CNV/{sample}/mosdepth.regions.bed.gz"
    params:
        prefix="results/04_CNV/{sample}/mosdepth"
    resources:
        runtime=config["resources"]["general"]["runtime"],
        mem_mb=config["resources"]["general"]["mem_mb"],
        cpus_per_task=config["resources"]["general"]["cpus"]
    container:
        "workflow/containers/mosdepth.sif"
    shell:
        """
        mkdir -p results/04_CNV/{wildcards.sample}

        mosdepth \
            --threads {resources.cpus_per_task} \
            --by {input.bed} \
            --no-per-base \
            {params.prefix} \
            {input.bam}
        """

rule normalize_region_depth:
    input:
        regions="results/04_CNV/{sample}/mosdepth.regions.bed.gz"
    output:
        tsv="results/04_CNV/{sample}/cnv_per_region.tsv"
    container:
        "workflow/containers/python.sif"
    resources:
        runtime=config["resources"]["general"]["runtime"],
        mem_mb=config["resources"]["general"]["mem_mb"],
        cpus_per_task=config["resources"]["general"]["cpus"]
    shell:
        """
        python workflow/scripts/normalize_cnv.py \
            --input {input.regions} \
            --output {output.tsv}
        """

rule merge_depth_cnv:
    input:
        tsvs=lambda _: expand(
            "results/04_CNV/{sample}/cnv_per_region.tsv",
            sample=get_passed_samples()
        )
    output:
        matrix="results/04_CNV/cnv_depth_matrix.tsv"
    container:
        "workflow/containers/python.sif"
    resources:
        runtime=config["resources"]["general"]["runtime"],
        mem_mb=config["resources"]["general"]["mem_mb"],
        cpus_per_task=config["resources"]["general"]["cpus"]
    shell:
        """
        sample_names=$(for f in {input.tsvs}; do basename $(dirname $f); done | tr '\n' ' ')

        python workflow/scripts/merge_cnv.py \
            --input-files {input.tsvs} \
            --sample-names $sample_names \
            --output {output.matrix}
        """

rule cnv_depth_to_presence:
    input:
        tsvs=lambda _: expand(
            "results/04_CNV/{sample}/cnv_per_region.tsv",
            sample=get_passed_samples()
        )
    output:
        matrix="results/cnvs_gene_presence.tsv"
    params:
        input_dir="results/04_CNV"
    resources:
        runtime=config["resources"]["general"]["runtime"],
        mem_mb=config["resources"]["general"]["mem_mb"],
        cpus_per_task=config["resources"]["general"]["cpus"]
    container:
        "workflow/containers/python.sif"
    shell:
        """
        sample_names=$(for f in {input.tsvs}; do basename $(dirname $f); done | tr '\n' ' ')

        python workflow/scripts/convert_cnv_to_matrix.py \
            --samples $sample_names  \
            --input-dir {params.input_dir} \
            --output {output.matrix}
        """
