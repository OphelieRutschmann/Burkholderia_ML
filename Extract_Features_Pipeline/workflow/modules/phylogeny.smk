# workflow/modules/phylogeny.smk
# Build phylogenetic tree from core SNP alignment
include: "../rules/common.smk"

rule build_phylogeny:
    input:
        aln="results/02_Variant_SNPs/core_SNP/clean.full.aln"
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
        iqtree2 -s {input.aln} \
               -m GTR+G \
               -nt AUTO \
               -pre {params.prefix}
        """
