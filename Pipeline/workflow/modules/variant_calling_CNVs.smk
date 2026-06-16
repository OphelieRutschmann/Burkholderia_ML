# workflow/modules/variant_calling_cnv_depth.smk

include: "../rules/common.smk"

rule prepare_regions:
    input:
        gff="databases/genomes/refgenome/ref_genome.gbff",
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

            name="unknown"
            id="unknown"

            n=split($9,a,";")

            for(i=1;i<=n;i++) {{

                if(a[i] ~ /^ID=/) {{
                    id=substr(a[i],4)
                }}

                if(a[i] ~ /^Name=/) {{
                    name=substr(a[i],6)
                }}
            }}

            feature = id"__"name"__"$3
            print $1,$4-1,$5,feature,".",$7
        }}' \
        > {output.bed}.tmp_features

        # Extract promoter regions for CDS and gene features only.
        # bedtools flank extends each feature upstream (strand-aware, -l bp, -r 0)
        # to approximate the promoter window defined in the config.
        grep -v "^#" {input.gff} | \
        awk '$3=="CDS" || $3=="gene"' | \
        awk 'BEGIN{{OFS="\t"}} {{

            name="unknown"
            id="unknown"

            n=split($9,a,";")

            for(i=1;i<=n;i++) {{

                if(a[i] ~ /^ID=/) {{
                    id=substr(a[i],4)
                }}

                if(a[i] ~ /^Name=/) {{
                    name=substr(a[i],6)
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
        cat \
            {output.bed}.tmp_features \
            {output.bed}.tmp_promoters | \
        bedtools sort -i stdin | \
        uniq > {output.bed}

        rm \
            {output.bed}.tmp_features \
            {output.bed}.tmp_promoters
        """

rule mosdepth_regions:
    input:
        bam="results/01_alignment/{sample}.bam",
        bai="results/01_alignment/{sample}.bam.bai",
        bed="databases/genomes/refgenome/regions.bed"
    output:
        regions="results/04_CNV_depth/{sample}/mosdepth.regions.bed.gz"
    params:
        prefix="results/04_CNV_depth/{sample}/mosdepth"
    resources:
        runtime=config["resources"]["general"]["runtime"],
        mem_mb=config["resources"]["general"]["mem_mb"],
        cpus_per_task=config["resources"]["general"]["cpus"]
    container:
        "workflow/containers/mosdepth.sif"
    shell:
        """
        mkdir -p results/04_CNV_depth/{wildcards.sample}

        mosdepth \
            --threads {resources.cpus_per_task} \
            --by {input.bed} \
            --no-per-base \
            {params.prefix} \
            {input.bam}
        """

rule normalize_region_depth:
    input:
        regions="results/04_CNV_depth/{sample}/mosdepth.regions.bed.gz"
    output:
        tsv="results/04_CNV_depth/{sample}/cnv_per_region.tsv"
    resources:  # FIX: missing resources block added for scheduler compatibility
        runtime=config["resources"]["general"]["runtime"],
        mem_mb=config["resources"]["general"]["mem_mb"],
        cpus_per_task=config["resources"]["general"]["cpus"]
    run:
        import numpy as np
        import pandas as pd

        df = pd.read_csv(
            input.regions,
            sep="\t",
            compression="gzip",
            header=None,
            names=["chrom", "start", "end", "region", "depth"]
        )

        # Use only CDS and gene features to estimate a stable per-chromosome
        # baseline coverage (median), excluding promoter regions which may have
        # systematically lower depth and would bias the normalisation.
        mask = df["region"].str.contains("__CDS") | df["region"].str.contains("__gene")

        medians = df.loc[mask].groupby("chrom")["depth"].median().to_dict()

        df["chrom_median"] = df["chrom"].map(medians)
        df["cnv_ratio"] = df["depth"] / df["chrom_median"]
        # Zeros replaced with NaN before log2 to avoid -inf in the feature matrix.
        df["log2cnv"] = np.log2(df["cnv_ratio"].replace(0, np.nan))
        df[["region", "chrom", "depth", "chrom_median", "cnv_ratio", "log2cnv"]].to_csv(
            output.tsv, sep="\t", index=False
        )


rule merge_depth_cnv:
    input:
        tsvs=lambda _: expand(
            "results/04_CNV_depth/{sample}/cnv_per_region.tsv",
            sample=get_passed_samples()
        )
    output:
        matrix="results/04_CNV_depth/cnv_depth_matrix.tsv"
    params:
        samples=lambda _: get_passed_samples()
    run:
        import pandas as pd

        dfs = []

        for sample, tsv in zip(params.samples, input.tsvs):
            df = pd.read_csv(tsv, sep="\t")
            df = df[["region", "log2cnv"]]
            df.columns = ["region", sample]
            df = df.set_index("region")
            dfs.append(df)

        # Concatenate per-sample Series along columns, then transpose so that
        # the final matrix is (samples × regions) — one row per sample, one
        # column per genomic feature, ready for ML ingestion.
        matrix = pd.concat(dfs, axis=1).T
        matrix.index.name = "sample"
        matrix.columns = [f"cnv_depth__{c}" for c in matrix.columns]
        matrix.to_csv(output.matrix, sep="\t")