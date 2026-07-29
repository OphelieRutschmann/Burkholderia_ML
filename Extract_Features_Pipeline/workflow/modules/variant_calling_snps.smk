# workflow/modules/variant_calling_snps.smk
include: "../rules/common.smk"

## This module performs SNP calling using snippy, followed by filtering and annotation of the vcf file using SNPEff ##

rule fix_bam_coordinates:
    input:
        bam="results/01_alignement/{sample}.aligned.sorted.bam",
        bai="results/01_alignement/{sample}.aligned.sorted.bam.bai",
        ref="databases/genomes/refgenome/ref_genome.fasta"
    output:
        fixed_bam="results/02_Variant_SNPs/{sample}/snps.fixed.bam",
        fixed_bai="results/02_Variant_SNPs/{sample}/snps.fixed.bam.bai"
    resources:
        runtime=config["resources"]["general"]["runtime"],
        mem_mb=config["resources"]["general"]["mem_mb"],
        cpus_per_task=config["resources"]["general"]["cpus"]
    container:
        "workflow/containers/samtools.sif"
    shell:
        """
        mkdir -p $(dirname {output.fixed_bam})

        # Reindex BAM with current reference to fix any coordinate mismatches
        samtools view -b -h -T {input.ref} {input.bam} > {output.fixed_bam}
        samtools index {output.fixed_bam}
        """

rule freebayes_call:
    input:
        bam="results/02_Variant_SNPs/{sample}/snps.fixed.bam",
        bai="results/02_Variant_SNPs/{sample}/snps.fixed.bam.bai",
        ref="databases/genomes/refgenome/ref_genome.fasta"
    output:
        vcf="results/02_Variant_SNPs/{sample}/snps.vcf"
    resources:
        runtime=config["resources"]["general"]["runtime"],
        mem_mb=config["resources"]["general"]["mem_mb"],
        cpus_per_task=config["resources"]["general"]["cpus"]
    container:
        "workflow/containers/freebayes.sif"
    shell:
        """
        mkdir -p $(dirname {output.vcf})

        freebayes -f {input.ref} \
            --min-coverage 10 \
            --min-alternate-fraction 0.9 \
            --min-mapping-quality 50 \
            --min-base-quality 50 \
            --ploidy 1 \
            {input.bam} > {output.vcf}
        """

rule filter_freebayes_vcf:
    input:
        vcf="results/02_Variant_SNPs/{sample}/snps.vcf"
    output:
        filt_vcf="results/02_Variant_SNPs/{sample}/snps.filt.vcf"
    resources:
        runtime=config["resources"]["general"]["runtime"],
        mem_mb=config["resources"]["general"]["mem_mb"],
        cpus_per_task=config["resources"]["general"]["cpus"]
    container:
        "workflow/containers/bcftools.sif"
    shell:
        """
        # Filter VCF: keep variants with minimum quality score
        bcftools filter -i 'QUAL>=50' {input.vcf} -o {output.filt_vcf}
        """

rule index_snippy_vcf:
    input:
        vcf="results/02_Variant_SNPs/{sample}/snps.filt.vcf"
    output:
        vcf_gz="results/02_Variant_SNPs/{sample}/snps.filt.vcf.gz",
        tbi="results/02_Variant_SNPs/{sample}/snps.filt.vcf.gz.tbi"
    resources:
        runtime=config["resources"]["general"]["runtime"], 
        mem_mb=config["resources"]["general"]["mem_mb"],    
        cpus_per_task=config["resources"]["general"]["cpus"] 
    container:
        "workflow/containers/bcftools.sif"
    shell:
        """
        # compress
        bcftools view -Oz -o {output.vcf_gz} {input.vcf}

        # index
        tabix -p vcf {output.vcf_gz}
        """

rule merge_filter_vcfs:
    input:
        lambda _: expand(
        "results/02_Variant_SNPs/{sample}/snps.filt.vcf.gz",
        sample=get_passed_samples()
        )
    output:
        merged_vcf="results/02_Variant_SNPs/merged.vcf.gz",
        merged_filtered_vcf="results/02_Variant_SNPs/merged.filt.vcf.gz"
    resources:
        runtime=config["resources"]["general"]["runtime"],
        mem_mb=config["resources"]["general"]["mem_mb"],
        cpus_per_task=config["resources"]["general"]["cpus"]
    container:
        "workflow/containers/bcftools.sif"
    shell:
        """
        bcftools merge {input} -Oz -o {output.merged_vcf} --missing-to-ref --force-single
        tabix -p vcf {output.merged_vcf}

        # remove variants where more than 10% of samples have a missing genotype.
        bcftools view -i 'F_MISSING < 0.1' {output.merged_vcf} -Oz -o {output.merged_filtered_vcf}

        # reindex to fix contig header warnings
        tabix -p vcf {output.merged_filtered_vcf}
        """

rule annotate_ref:
    input:
        ref="databases/genomes/refgenome/ref_genome.fasta",
        flag="databases/bakta_db/.download_complete"
    output:
        ref_annotation="databases/genomes/refgenome/ref_genome.gbff"
    params:
        outdir="databases/genomes/refgenome",
        db_path=config["databases_bakta"]["bakta_path"]
    resources:
        runtime=config["resources"]["general"]["runtime"], 
        mem_mb=config["resources"]["general"]["mem_mb"],    
        cpus_per_task=config["resources"]["general"]["cpus"] 
    container:
        "workflow/containers/bakta.sif"
    shell:
        """
        echo "annotating reference strain {input.ref}"

        if [ -n "{params.db_path}" ] && [ -d "{params.db_path}" ]; then
            BAKTA_DB="{params.db_path}"
        else
            BAKTA_DB="databases/bakta_db/db-light"
        fi

        bakta \
            --db $BAKTA_DB \
            --output {params.outdir} \
            --force \
            --prefix ref_genome \
            --keep-contig-headers \
            {input.ref}
        """

rule build_snpeff_db:
    input:
        ref="databases/genomes/refgenome/ref_genome.fasta",
        gbff="databases/genomes/refgenome/ref_genome.gbff"
    output:
        db_done="databases/snpeff_db/.db_built"
    params:
        outdir="databases/snpeff_db"
    resources:
        runtime=config["resources"]["general"]["runtime"],
        mem_mb=config["resources"]["general"]["mem_mb"],
        cpus_per_task=1
    container:
        "workflow/containers/snpEff.sif"
    shell:
        """
        mkdir -p {params.outdir}/data/ref
        
        # Symlink files to snpEff expected structure
        ln -sf $(realpath {input.gbff}) {params.outdir}/data/ref/genes.gbk
        ln -sf $(realpath {input.ref}) {params.outdir}/data/ref/sequences.fa
        
        # Create snpEff config
        echo "data.dir = $(realpath {params.outdir})/data" > {params.outdir}/snpEff.config
        echo "ref.genome : ref" >> {params.outdir}/snpEff.config
        
        # Build database
        
        snpEff build -genbank -v \
            -c {params.outdir}/snpEff.config \
            -noCheckCds \
            -noCheckProtein \
            ref
        
        touch {output.db_done}
        """

rule snpeff:
    input:
        vcf="results/02_Variant_SNPs/merged.filt.vcf.gz",
        db_done="databases/snpeff_db/.db_built"
    output:
        vcf="results/02_Variant_SNPs/snps.ann.vcf",
        stats="results/02_Variant_SNPs/snpeff_stats.html"
    params:
        db_dir="databases/snpeff_db"
    resources:
        runtime=config["resources"]["general"]["runtime"],
        mem_mb=config["resources"]["general"]["mem_mb"],
        cpus_per_task=config["resources"]["general"]["cpus"]
    container:
        "workflow/containers/snpEff.sif"
    shell:
        """
        snpEff ann -v \
            -c {params.db_dir}/snpEff.config \
            -stats {output.stats} \
            ref \
            {input.vcf} > {output.vcf}

        # Verify output VCF structure (fixes contig header warnings)
        echo "SNPEff annotation complete"
        """

# ============================================================================
# Core SNP extraction and feature matrix generation
# ============================================================================

rule core_snp_alignment:
    input:
        ref="databases/genomes/refgenome/ref_genome.fasta",
        snippy_vcfs=expand(
            "results/02_Variant_SNPs/{sample}/snps.filt.vcf",
            sample=get_passed_samples()
        )
    output:
        core_aln="results/02_Variant_SNPs/core_SNP/core.full.aln",
        clean_aln="results/02_Variant_SNPs/core_SNP/clean.full.aln",
        snp_pos="results/02_Variant_SNPs/core_SNP/snp_positions.txt"
    params:
        snippy_dir="results/02_Variant_SNPs",
        outdir="results/02_Variant_SNPs/core_SNP",
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

rule core_snp_feature_matrix:
    input:
        aln="results/02_Variant_SNPs/core_SNP/clean.full.aln"
    output:
        matrix="results/snps_core.tsv"
    container:
        "workflow/containers/biopython.sif"
    shell:
        """
        python workflow/scripts/convert_snp_to_matrix.py \
            --alignment {input.aln} \
            --output {output.matrix}
        """