# workflow/modules/phylogeny.smk
# Build phylogenetic tree from core SNP alignment

include: "../rules/common.smk"

rule remove_recombination:
    input:
        aln="results/02_Variant_SNPs/core_SNP/clean.full.aln"
    output:
        filtered_aln="results/04_Phylogeny/core_SNP/clean.filtered_polymorphic_sites.fasta",
        recomb_gff="results/04_Phylogeny/core_SNP/clean.recombination_predictions.gff"
    params:
        outdir="results/04_Phylogeny/core_SNP",
        prefix="clean"
    threads:
        config["resources"]["general"]["cpus"]
    resources:
        runtime=config["resources"]["general"]["runtime"],
        mem_mb=config["resources"]["general"]["mem_mb"],
        cpus_per_task=config["resources"]["general"]["cpus"]
    container:
        "workflow/containers/gubbins.sif"
    shell:
        """
        mkdir -p {params.outdir}

        run_gubbins.py --threads {threads} \
                       --prefix {params.outdir}/{params.prefix} \
                       --verbose \
                       --filter-percentage 25 \
                       {input.aln}
        """

rule build_phylogeny:
    input:
        aln="results/04_Phylogeny/core_SNP/clean.filtered_polymorphic_sites.fasta"
    output:
        tree="results/04_Phylogeny/iqtree/iqtree.treefile",
        model="results/04_Phylogeny/iqtree/iqtree.model.gz"
    params:
        prefix="results/04_Phylogeny/iqtree/iqtree"
    resources:
        runtime=config["resources"]["general"]["runtime"],
        mem_mb=config["resources"]["general"]["mem_mb"],
        cpus_per_task=config["resources"]["general"]["cpus"]
    container:
        "workflow/containers/iqtree.sif"
    shell:
        """
        mkdir -p $(dirname {params.prefix})

        iqtree -s {input.aln} \
               -m TEST \
               -bb 1000 \
               -nt AUTO \
               -pre {params.prefix}
        """