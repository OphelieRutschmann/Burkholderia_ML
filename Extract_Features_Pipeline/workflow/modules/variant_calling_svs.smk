# workflow/modules/variant_calling_svs.smk
include: "../rules/common.smk"

## This module detects structural variants using Pindel. 
## Before using Pindel, PCR duplicates are marked using Picard.
## Per-sample calls are merged across samples with SURVIVOR (breakpoint distance
## tolerance), and SV loci are checked for read coverage in every sample with
## mosdepth so that SVs lying in accessory regions are excluded from the matrix.

# Optional settings (add an `sv:` section to config/config.yaml to override)
SV_CFG = config.get("sv", {})
SV_MAX_DIST = SV_CFG.get("max_dist", 100)      # max breakpoint distance (bp) to merge calls across samples
SV_MIN_SIZE = SV_CFG.get("min_size", 50)       # minimum SV length (bp) kept by SURVIVOR
SV_FLANK = SV_CFG.get("flank", 500)            # flank size (bp) checked on each side of an SV
SV_MIN_DEPTH = SV_CFG.get("min_depth", 10)     # depth a flank base needs to count as covered
SV_MIN_FRAC = SV_CFG.get("min_frac", 0.9)      # fraction of flank bases that must reach SV_MIN_DEPTH

# Use alignement to reference generated in the QC_Illumina.smk module (for coverage calculation)
rule mark_duplicates:
    input:
        bam="results/01_alignement/{sample}.aligned.sorted.bam"
    output:
        bam="results/01_alignement/{sample}.dedup.bam",
        bai="results/01_alignement/{sample}.dedup.bai",
        metrics="results/01_alignement/{sample}.dedup.metrics.txt"
    resources:
        runtime=config["resources"]["general"]["runtime"],
        mem_mb=config["resources"]["general"]["mem_mb"],
        cpus_per_task=config["resources"]["general"]["cpus"]
    container:
        "workflow/containers/picard.sif"
    shell:
            """
            java -jar /usr/picard/picard.jar MarkDuplicates \
                I={input.bam} \
                O={output.bam} \
                M={output.metrics} \
                REMOVE_DUPLICATES=false \
                CREATE_INDEX=true \
                VALIDATION_STRINGENCY=SILENT
            """

rule pindel_config:
    input:
        bam="results/01_alignement/{sample}.dedup.bam"
    output:
        config="results/03_SV/{sample}/pindel_config.txt",
        metrics="results/03_SV/{sample}/insert_metrics.txt"
    params:
        outdir="results/03_SV/{sample}"
    resources:
        runtime=config["resources"]["general"]["runtime"],
        mem_mb=config["resources"]["general"]["mem_mb"],
        cpus_per_task=config["resources"]["general"]["cpus"]
    container:
        "workflow/containers/picard.sif"
    shell:
        """
        mkdir -p {params.outdir}
        
        # Compute insert size metrics
        java -jar /usr/picard/picard.jar CollectInsertSizeMetrics \
            I={input.bam} \
            O={output.metrics} \
            H={params.outdir}/insert_histogram.pdf

        # Extract the mean insert size (what Pindel expects): locate the
        # MEAN_INSERT_SIZE column in the header line, then read it from the first data line.
        INSERT_SIZE=$(awk 'BEGIN{{FS="\t"}} /^#/{{next}} /^MEDIAN_INSERT_SIZE/{{for(i=1;i<=NF;i++) if($i=="MEAN_INSERT_SIZE") c=i; next}} c && NF>0{{print int($c); exit}}' {output.metrics})
        if [ -z "$INSERT_SIZE" ] || [ "$INSERT_SIZE" -le 0 ]; then
            echo "Could not read MEAN_INSERT_SIZE from {output.metrics}" >&2
            exit 1
        fi

        echo "insert size: $INSERT_SIZE"

        echo "{input.bam} $INSERT_SIZE {wildcards.sample}" > {output.config}
        """
        
rule fasta_index:
    input:
        ref="databases/genomes/refgenome/ref_genome.fasta"
    output:
        fai="databases/genomes/refgenome/ref_genome.fasta.fai"
    resources:
        runtime=config["resources"]["general"]["runtime"],
        mem_mb=config["resources"]["general"]["mem_mb"],
        cpus_per_task=config["resources"]["general"]["cpus"]
    container:
        "workflow/containers/samtools.sif"
    shell:
        """
        samtools faidx {input.ref}
        """

rule pindel:
    input:
        config="results/03_SV/{sample}/pindel_config.txt",
        ref="databases/genomes/refgenome/ref_genome.fasta",
        fai="databases/genomes/refgenome/ref_genome.fasta.fai"
    output:
        td="results/03_SV/{sample}/pindel_TD",
        inv="results/03_SV/{sample}/pindel_INV",
        li="results/03_SV/{sample}/pindel_LI",
        d="results/03_SV/{sample}/pindel_D",
        flag=touch("results/03_SV/{sample}/.pindel_done")
    params:
        prefix="results/03_SV/{sample}/pindel"
    resources:
        runtime=config["resources"]["general"]["runtime"],
        mem_mb=config["resources"]["general"]["mem_mb"],
        cpus_per_task=config["resources"]["general"]["cpus"]
    container:
        "workflow/containers/pindel.sif"
    shell:
        """
        pindel \
            -f {input.ref} \
            -i {input.config} \
            -o {params.prefix} \
            -T {resources.cpus_per_task}
        """

rule pindel2vcf:
    input:
        ref="databases/genomes/refgenome/ref_genome.fasta",
        flag="results/03_SV/{sample}/.pindel_done",
        pindel="results/03_SV/{sample}/pindel_{svtype}"
    output:
        vcf="results/03_SV/{sample}/sv_{svtype}.vcf"
    params:
        refname="ref",
        refdate="20240101"
    resources:
        runtime=config["resources"]["general"]["runtime"],
        mem_mb=config["resources"]["general"]["mem_mb"],
        cpus_per_task=config["resources"]["general"]["cpus"]
    container:
        "workflow/containers/pindel.sif"
    shell:
        """
        pindel2vcf \
            -p {input.pindel} \
            -r {input.ref} \
            -R {params.refname} \
            -d {params.refdate} \
            -v {output.vcf}
        """

rule index_sv_vcf:
    input:
        vcf="results/03_SV/{sample}/sv_{svtype}.vcf",
        fai="databases/genomes/refgenome/ref_genome.fasta.fai"
    output:
        vcf_gz="results/03_SV/{sample}/sv_{svtype}.vcf.gz",
        tbi="results/03_SV/{sample}/sv_{svtype}.vcf.gz.tbi"
    resources:
        runtime=config["resources"]["general"]["runtime"],
        mem_mb=config["resources"]["general"]["mem_mb"],
        cpus_per_task=config["resources"]["general"]["cpus"]
    container:
        "workflow/containers/bcftools.sif"
    shell:
        """
        # pindel2vcf writes no ##contig header lines, which bcftools needs: add them from the
        # reference index. Its output is also not guaranteed to be coordinate-sorted.
        bcftools reheader --fai {input.fai} {input.vcf} | bcftools sort -Oz -o {output.vcf_gz} -
        bcftools index -t {output.vcf_gz}
        """

# Combine the 4 SV types of one sample into a single per-sample VCF.
# (Merging them directly with bcftools merge would create 4 columns per sample.)
# Combine the 4 SV types of one sample into a single, sorted, uncompressed VCF
# (SURVIVOR does not read bgzipped input).
# Pindel reports deletions with inserted sequence as SVTYPE=RPL, a type SURVIVOR
# does not recognise, so they are relabelled as DEL.
rule concat_sample_svs:
    input:
        vcfs=expand("results/03_SV/{{sample}}/sv_{svtype}.vcf.gz", svtype=["TD", "INV", "LI", "D"]),
        tbis=expand("results/03_SV/{{sample}}/sv_{svtype}.vcf.gz.tbi", svtype=["TD", "INV", "LI", "D"])
    output:
        vcf="results/03_SV/{sample}/sv_all.vcf"
    resources:
        runtime=config["resources"]["general"]["runtime"],
        mem_mb=config["resources"]["general"]["mem_mb"],
        cpus_per_task=config["resources"]["general"]["cpus"]
    container:
        "workflow/containers/bcftools.sif"
    shell:
        """
        bcftools concat -a -D {input.vcfs} -Ov | \
            sed 's/SVTYPE=RPL/SVTYPE=DEL/' | \
            bcftools sort -Ov -o {output.vcf}
        """

# Merge SVs across samples with SURVIVOR: calls of the same type whose breakpoints
# are within SV_MAX_DIST bp are treated as the same event.
# SURVIVOR arguments: list, max distance, min support (1 sample), same type required (1),
# strand not required (0), no distance estimation (0), min SV size.
# The INFO/SUPP_VEC field gives presence (1) per sample, in the order of the input list.
rule merge_sv_vcfs:
    input:
        vcfs=lambda _: expand("results/03_SV/{sample}/sv_all.vcf", sample=get_passed_samples())
    output:
        raw=temp("results/03_SV/merged_svs.survivor.vcf"),
        samples="results/03_SV/merged_svs.samples.txt"
    params:
        vcf_list="results/03_SV/survivor_input.txt",
        samples=lambda _: get_passed_samples(),
        max_dist=SV_MAX_DIST,
        min_size=SV_MIN_SIZE
    resources:
        runtime=config["resources"]["general"]["runtime"],
        mem_mb=config["resources"]["general"]["mem_mb"],
        cpus_per_task=config["resources"]["general"]["cpus"]
    container:
        "workflow/containers/survivor.sif"
    shell:
        """
        printf '%s\n' {input.vcfs} > {params.vcf_list}
        printf '%s\n' {params.samples} > {output.samples}

        SURVIVOR merge {params.vcf_list} {params.max_dist} 1 1 0 0 {params.min_size} {output.raw}
        """

# Sort, give every merged SV a unique ID (sv__CHROM_POS_END_TYPE) and compress.
rule index_merged_svs:
    input:
        raw="results/03_SV/merged_svs.survivor.vcf",
        samples="results/03_SV/merged_svs.samples.txt",
        fai="databases/genomes/refgenome/ref_genome.fasta.fai"
    output:
        merged="results/03_SV/merged_svs.vcf.gz",
        tbi="results/03_SV/merged_svs.vcf.gz.tbi"
    resources:
        runtime=config["resources"]["general"]["runtime"],
        mem_mb=config["resources"]["general"]["mem_mb"],
        cpus_per_task=config["resources"]["general"]["cpus"]
    container:
        "workflow/containers/bcftools.sif"
    shell:
        """
        bcftools reheader --fai {input.fai} {input.raw} | bcftools sort -Ov - | \
        awk 'BEGIN{{FS=OFS="\t"}} /^#/{{print; next}} {{
            end=$2; t="NA"
            n=split($8, a, ";")
            for (i=1; i<=n; i++) {{
                if (a[i] ~ /^END=/) end=substr(a[i], 5)
                if (a[i] ~ /^SVTYPE=/) t=substr(a[i], 8)
            }}
            id="sv__"$1"_"$2"_"end"_"t
            if (id in seen) {{ seen[id]++; id=id"_"seen[id] }} else {{ seen[id]=1 }}
            $3=id
            print
        }}' | \
        bcftools view -Oz -o {output.merged}
        bcftools index -t {output.merged}
        """

# BED of the regions flanking each merged SV (left of POS, right of END), clipped to
# the contig. Flanks lie outside the SV itself, so they stay covered in a sample that
# carries a deletion, but are uncovered when the whole locus is absent from that strain.
rule sv_flanks:
    input:
        vcf="results/03_SV/merged_svs.vcf.gz",
        fai="databases/genomes/refgenome/ref_genome.fasta.fai"
    output:
        bed="results/03_SV/sv_flanks.bed"
    params:
        flank=SV_FLANK
    resources:
        runtime=config["resources"]["general"]["runtime"],
        mem_mb=config["resources"]["general"]["mem_mb"],
        cpus_per_task=config["resources"]["general"]["cpus"]
    container:
        "workflow/containers/bcftools.sif"
    shell:
        """
        bcftools query -f '%CHROM\t%POS\t%INFO/END\t%ID\n' {input.vcf} | \
        awk -v F={params.flank} 'BEGIN{{OFS="\t"}} NR==FNR{{len[$1]=$2; next}} {{
            s=($2 < $3) ? $2 : $3
            e=($2 < $3) ? $3 : $2
            ls=s-1-F; if (ls < 0) ls=0
            if (ls < s-1) print $1, ls, s-1, $4"__L"
            re=e+F; if (re > len[$1]) re=len[$1]
            if (e < re) print $1, e, re, $4"__R"
        }}' {input.fai} - | sort -k1,1 -k2,2n > {output.bed}
        """

# Per-sample coverage of the SV flanks: number of bases with depth >= SV_MIN_DEPTH.
rule mosdepth_sv_flanks:
    input:
        bam="results/01_alignement/{sample}.dedup.bam",
        bai="results/01_alignement/{sample}.dedup.bai",
        bed="results/03_SV/sv_flanks.bed"
    output:
        thresholds="results/03_SV/{sample}/flanks.thresholds.bed.gz"
    params:
        prefix="results/03_SV/{sample}/flanks",
        min_depth=SV_MIN_DEPTH
    resources:
        runtime=config["resources"]["general"]["runtime"],
        mem_mb=config["resources"]["general"]["mem_mb"],
        cpus_per_task=config["resources"]["general"]["cpus"]
    container:
        "workflow/containers/mosdepth.sif"
    shell:
        """
        mosdepth \
            --threads {resources.cpus_per_task} \
            --by {input.bed} \
            --thresholds {params.min_depth} \
            --no-per-base \
            {params.prefix} \
            {input.bam}
        """

rule extract_core_svs:
    input:
        vcf="results/03_SV/merged_svs.vcf.gz",
        samples="results/03_SV/merged_svs.samples.txt",
        thresholds=lambda _: expand("results/03_SV/{sample}/flanks.thresholds.bed.gz", sample=get_passed_samples())
    output:
        matrix="results/svs_core.tsv",
        metadata="results/svs_core_meta.tsv"
    params:
        coverage=lambda _: [f"{s}=results/03_SV/{s}/flanks.thresholds.bed.gz" for s in get_passed_samples()],
        min_frac=SV_MIN_FRAC
    resources:
        runtime=config["resources"]["general"]["runtime"],
        mem_mb=config["resources"]["general"]["mem_mb"],
        cpus_per_task=config["resources"]["general"]["cpus"]
    container:
        "workflow/containers/biopython.sif"
    shell:
        """
        python workflow/scripts/convert_sv_to_matrix.py \
            --vcf {input.vcf} \
            --samples {input.samples} \
            --coverage {params.coverage} \
            --min-frac {params.min_frac} \
            --output {output.matrix} \
            --metadata {output.metadata}
        """
