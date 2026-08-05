# workflow/modules/variant_calling_snps.smk
include: "../rules/common.smk"

## This module performs SNP calling using snippy, followed by filtering and annotation of the vcf file using SNPEff ##

rule variant_snippy:
    input:
        r1="results/00_QC/fastp/trimmed_reads/{sample}_trimmed_R1.fastq",
        r2="results/00_QC/fastp/trimmed_reads/{sample}_trimmed_R2.fastq",
        ref="databases/genomes/refgenome/ref_genome.fasta"
    output:
        vcf = "results/02_Variant_SNPs/{sample}/snps.filt.vcf",
        alignement = "results/02_Variant_SNPs/{sample}/snps.aligned.fa"
    params:
        outdir="results/02_Variant_SNPs/{sample}",
        read_type=config["read_type"]
    resources:
        runtime=config["resources"]["general"]["runtime"],
        mem_mb=config["resources"]["general"]["mem_mb"],
        cpus_per_task=config["resources"]["general"]["cpus"]
    container:
        "workflow/containers/snippy.sif"
    shell:
        """       
        rm -rf {params.outdir}
        mkdir -p {params.outdir}/tmp

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

        # Rename output to expected name
        mv "{params.outdir}/snps.vcf" "{params.outdir}/snps.filt.vcf"

        # Clean up temporary directory and other unnecessary files
        rm -rf {params.outdir}/tmp
        rm -rf {params.outdir}/reference
        rm -f "{params.outdir}/snps.html"
        rm -f "{params.outdir}/snps.bam"
        rm -f "{params.outdir}/snps.bam.bai"
        rm -f "{params.outdir}/snps.consensus.fa"
        rm -f "{params.outdir}/snps.consensus.subs.fa"
        rm -f "{params.outdir}/snps.bed"
        rm -f "{params.outdir}/snps.log"
        rm -f "{params.outdir}/snps.txt"
        rm -f "{params.outdir}/snps.tab"
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

rule core_snp_alignment:
    input:
        ref="databases/genomes/refgenome/ref_genome.fasta",
        checkpoint="results/00_QC/passed_samples.csv",
        snippy_vcfs=lambda _: expand(
            "results/02_Variant_SNPs/{sample}/snps.filt.vcf",
            sample=get_passed_samples()
        ),
        snippy_alignements=lambda _: expand(
            "results/02_Variant_SNPs/{sample}/snps.aligned.fa",
            sample=get_passed_samples()
        ),
        snippy_dirs=lambda _: expand(
            "results/02_Variant_SNPs/{sample}",
            sample=get_passed_samples()
        )
    output:
        core_aln="results/02_Variant_SNPs/core_SNP/core.full.aln",
        clean_aln="results/02_Variant_SNPs/core_SNP/clean.full.aln",
        snp_pos="results/02_Variant_SNPs/core_SNP/snp_positions.txt"
    params:
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

	for d in {input.snippy_dirs}; do
	    cp "$d/snps.filt.vcf" "$d/snps.vcf"
	done


        # Create core alignment from all snippy VCFs
        snippy-core \
            --ref {input.ref} \
            --prefix {params.outdir}/{params.prefix} \
            {input.snippy_dirs}

        # Remove reference sequence from alignment
	ref_header=$(grep "^>" {input.ref} | head -1 | sed 's/^>//')
	snippy-clean_full_aln {params.outdir}/{params.prefix}.full.aln \
	    | awk -v ref=">$ref_header" '
		/^>/ {{ skip = ($0 == ref) }}
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
        matrix="results/02_Variant_SNPs/core_SNP/core.tsv"
    resources:
        runtime=config["resources"]["general"]["runtime"],
        mem_mb=config["resources"]["general"]["mem_mb"],
        cpus_per_task=config["resources"]["general"]["cpus"]
    container:
        "workflow/containers/biopython.sif"
    shell:
        """
        python workflow/scripts/convert_snp_to_matrix.py \
            --alignment {input.aln} \
            --output {output.matrix}
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

rule core_snps_to_vcf:
    input:
        core_tab="results/02_Variant_SNPs/core_SNP/core.tsv",
        ref="databases/genomes/refgenome/ref_genome.fasta"
    output:
        vcf="results/02_Variant_SNPs/core_SNP/core_snps.vcf"
    resources:
        runtime=config["resources"]["general"]["runtime"],
        mem_mb=config["resources"]["general"]["mem_mb"],
        cpus_per_task=config["resources"]["general"]["cpus"]
    container:
        "workflow/containers/biopython.sif"
    shell:
        """
        python workflow/scripts/core_snps_to_vcf.py \
            {input.core_tab} \
            {input.ref} \
            {output.vcf}
        """

rule snpeff_core_snps:
    input:
        vcf="results/02_Variant_SNPs/core_SNP/core_snps.vcf",
        db_done="databases/snpeff_db/.db_built"
    output:
        vcf="results/02_Variant_SNPs/core_SNP/core_snps.ann.vcf",
        stats="results/02_Variant_SNPs/core_SNP/snpeff_stats.html"
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
