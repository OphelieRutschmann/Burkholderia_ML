# workflow/modules/variant_calling_snps.smk
include: "../rules/common.smk"

## This module performs SNP calling using snippy, followed by filtering and annotation of the vcf file using SNPEff ##

rule variant_snippy: 
    input:
        r1="results/00_QC/fastp/trimmed_reads/{sample}_trimmed_R1.fastq",
        r2="results/00_QC/fastp/trimmed_reads/{sample}_trimmed_R2.fastq",
        ref="databases/genomes/refgenome/ref_genome.fasta"
    output:
        "results/Variant_SNPs/{sample}/snps.filt.vcf"
    params:
        outdir="results/Variant_SNPs/{sample}",
        read_type=config["read_type"]
    resources:
        runtime=config["resources"]["general"]["runtime"], 
        mem_mb=config["resources"]["general"]["mem_mb"],    
        cpus_per_task=config["resources"]["general"]["cpus"] 
    container:
        "workflow/containers/snippy.sif"
    shell:
        """
        mkdir -p {params.outdir}/tmp

        #perform snippy
        snippy --cpus {resources.cpus_per_task} \
            --outdir {params.outdir} \
            --reference {input.ref} \
            --R1 {input.r1} \
            --R2 {input.r2} \
            --mincov 10 \
            --minfrac 0.9 \
            --minqual 50 \
            --force \
            --cleanup

        # Clean up temporary directory and other unecessary files
        rm -rf {params.outdir}/tmp
        rm -rf {params.outdir}/reference
        rm -f "{params.outdir}/snps.html"
        rm -f "{params.outdir}/snps.bam"
        rm -f "{params.outdir}/snps.bam.bai"
        rm -f "{params.outdir}/snps.consensus.fa"
        rm -f "{params.outdir}/snps.consensus.subs.fa"
        rm -f "{params.outdir}/snps.aligned.fa"
        rm -f "{params.outdir}/snps.bed"
        rm -f "{params.outdir}/snps.log"
        rm -f "{params.outdir}/snps.txt"
        rm -f "{params.outdir}/snps.tab"
        """

rule index_snippy_vcf:
    input:
        vcf="results/Variant_SNPs/{sample}/snps.filt.vcf"
    output:
        vcf_gz="results/Variant_SNPs/{sample}/snps.filt.vcf.gz",
        tbi="results/Variant_SNPs/{sample}/snps.filt.vcf.gz.tbi"
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
        "results/Variant_SNPs/{sample}/snps.filt.vcf.gz",
        sample=get_passed_samples()
        )
    output:
        merged_vcf="results/Variant_SNPs/merged.vcf.gz",
        merged_filtered_vcf="results/Variant_SNPs/merged.filt.vcf.gz"
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

        # Filter by allele frequency and missingness:
        # Remove variants (alleles) present in less than 1% of the samples, and remove sites where > 10% of the data is missing
        bcftools view -i 'MAF > 0.01 && F_MISSING < 0.1' {output.merged_vcf} -Oz -o {output.merged_filtered_vcf}

        # Index output
        tabix -p vcf {output.merged_filtered_vcf}
        """

rule annotate_ref:
    input:
        ref="databases/genomes/refgenome/ref_genome.fasta",
        flag="databases/bakta_db/.download_complete"
    output:
        ref_annotation="databases/genomes/refgenome/ref_genome.gff3"
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
            --prefix reference \
            {input.ref}
        """

rule build_snpeff_db:
    input:
        ref="databases/genomes/refgenome/ref_genome.fasta",
        gff="databases/genomes/refgenome/ref_genome.gff3"
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
        mkdir -p {params.outdir}
        
        # Copy files to snpEff database structure
        cp {input.gff} {params.outdir}/genes.gff
        
        # Create snpEff config
        echo "ref.genome : ref" > {params.outdir}/snpEff.config
        
        # Build database using the Bakta annotation
        snpEff build -gff3 -v -c {params.outdir}/snpEff.config ref \
        -noCheckCds -noCheckProtein 
               
        touch {output.db_done}
        """

rule snpeff:
    input:
        vcf="results/Variant_SNPs/merged.vcf.gz",
        db_done="databases/snpeff_db/.db_built"
    output:
        vcf="results/Variant_SNPs/snps.ann.vcf",
        stats="results/Variant_SNPs/snpeff_stats.html"
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
        """