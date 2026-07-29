# workflow/modules/phylogeny.smk
# Build core-SNP phylogeny from existing snippy outputs

include: "../rules/common.smk"

rule core_snp_alignment:
    input:
        ref="databases/genomes/refgenome/ref_genome.fasta",
        snippy_vcfs=expand(
            "results/02_Variant_SNPs/{sample}/snps.filt.vcf",
            sample=get_passed_samples()
        )
    output:
        core_aln="results/04_Phylogeny/core_SNP/core.full.aln",
        clean_aln="results/04_Phylogeny/core_SNP/clean.full.aln",
        snp_pos="results/04_Phylogeny/core_SNP/snp_positions.txt"
    params:
        snippy_dir="results/02_Variant_SNPs",
        outdir="results/04_Phylogeny/core_SNP",
        prefix="core"
    resources:
        runtime=config["resources"]["general"]["runtime"],
        mem_mb=config["resources"]["general"]["mem_mb"],
        cpus_per_task=config["resources"]["general"]["cpus"]
    container:
        "workflow/containers/snippy.sif"
    shell:
        """
        mkdir -p {params.outdir}

        # Create core alignment from all snippy VCFs
        snippy-core \
            --ref {input.ref} \
            --prefix {params.outdir}/{params.prefix} \
            {params.snippy_dir}/*/

        # Remove reference sequence from alignment
        ref_header=$(grep "^>" {input.ref} | head -1 | sed 's/^>//')
        snippy-clean_full_aln {params.outdir}/{params.prefix}.full.aln \
            | awk -v ref="$ref_header" '
                /^>/ {{ skip = ($0 == ">" ref); next }}
                !skip
            ' > {output.clean_aln}

        # Extract SNP positions
        cut -d$'\t' -f2 "{params.outdir}/core.tab" | tail -n +2 > "{output.snp_pos}"

        echo "Core alignment: {output.clean_aln}"
        echo "Number of SNP positions: $(wc -l < {output.snp_pos})"
        """

rule remove_recombination:
    input:
        aln="results/04_Phylogeny/core_SNP/clean.full.aln"
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