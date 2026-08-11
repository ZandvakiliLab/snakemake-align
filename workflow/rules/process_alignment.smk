rule samtools_collate:
    input:
        get_bam,
    output:
        temp("results/processed_alignment/collate/{sample}.bam"),
    log:
        "results/processed_alignment/collate/{sample}.log",
    threads: 2
    params:
        extra="-f",
    wrapper:
        "v9.15.0/bio/samtools/collate"


rule samtools_fixmate:
    # Only necessary if using samtools markdup for deduplication, as it requires correct mate information
    input:
        rules.samtools_collate.output,
    output:
        temp("results/processed_alignment/fixmate/{sample}.bam"),
    log:
        "results/processed_alignment/fixmate/{sample}.log",
    threads: 1
    params:
        extra="-m",
    message:
        "Fixing mate information in {wildcards.sample}"
    wrapper:
        "v9.15.0/bio/samtools/fixmate/"


rule samtools_sort:
    input:
        get_bam_2,
    output:
        temp("results/processed_alignment/sort/{sample}.bam"),
    log:
        "results/processed_alignment/sort/{sample}.log",
    threads: 2
    params:
        extra=config["mapping_postprocessing"]["samtools_sort"]["extra"],
    message:
        "re-sort reads after mapping regardless if mapper did"
    wrapper:
        "v9.4.1/bio/samtools/sort"


rule samtools_index:
    input:
        rules.samtools_sort.output,
    output:
        temp("results/processed_alignment/sort/{sample}.bam.bai"),
    log:
        "results/processed_alignment/sort/{sample}_index.log",
    threads: 2
    params:
        extra=config["mapping_postprocessing"]["samtools_index"]["extra"],
    message:
        "index reads"
    wrapper:
        "v9.4.1/bio/samtools/index"


########################
# DEDUPLICATION RULES
########################


rule umi_tools_dedup:
    input:
        bam=rules.samtools_sort.output,
        bai=rules.samtools_index.output,
    output:
        temp("results/processed_alignment/dedup/umi_tools/{sample}.bam"),
    log:
        "results/processed_alignment/dedup/umi_tools/{sample}.log",
    container:
        "docker://quay.io/biocontainers/umi_tools:1.1.6--py310h1fe012e_0"
    threads: 5
    params:
        extra=config["mapping_postprocessing"]["deduplication"]["umi_tools"]["extra"],
        paired="--paired" if is_paired_end() else "",
    message:
        "deduplicate reads using umi_tools"
    shell:
        """
        umi_tools dedup \
            -I {input} \
            -S {output} \
            --log={log} \
            {params.paired} \
            {params.extra}
        """


rule samtools_markdup:
    input:
        aln=rules.samtools_sort.output,
    output:
        bam=temp("results/processed_alignment/dedup/samtools/{sample}.bam"),
        metrics="results/processed_alignment/dedup/samtools/{sample}.txt",
    log:
        "results/processed_alignment/dedup/samtools/{sample}.log",
    threads: 2
    params:
        extra=config["mapping_postprocessing"]["deduplication"]["samtools"]["extra"],
    wrapper:
        "v9.15.0/bio/samtools/markdup"


rule add_RG:
    input:
        rules.samtools_sort.output,
    output:
        temp("results/processed_alignment/dedup/picard/{sample}_addRG.bam"),
    log:
        "results/processed_alignment/dedup/picard/{sample}_addRG.log",
    resources:
        mem_mb=1024,
    params:
        extra="--RGLB lib1 --RGPL illumina --RGPU {sample} --RGSM {sample}",
    wrapper:
        "v9.15.0/bio/picard/addorreplacereadgroups"


rule markduplicates_bam:
    input:
        bams=rules.add_RG.output,
    output:
        bam=temp("results/processed_alignment/dedup/picard/{sample}.bam"),
        metrics="results/processed_alignment/dedup/picard/{sample}.metrics.txt",
    log:
        "results/processed_alignment/dedup/picard/{sample}.log",
    resources:
        mem_mb=1024,
    params:
        extra="--REMOVE_DUPLICATES true",
    wrapper:
        "v9.15.0/bio/picard/markduplicates"


rule samtools_index_dedup:
    input:
        get_processed_bam,
    output:
        temp("results/processed_alignment/dedup/{tool}/{sample}.bam.bai"),
    log:
        "results/processed_alignment/dedup/{tool}/{sample}_index.log",
    threads: 2
    params:
        extra=config["mapping_postprocessing"]["samtools_index"]["extra"],
    message:
        "index reads"
    wrapper:
        "v9.4.1/bio/samtools/index"


########################
# BAM TO CRAM RULES
########################


rule bam_to_cram:
    input:
        bam=get_processed_bam,
        fa=rules.get_genome.output.fasta,
    output:
        "results/processed_alignment/cram/{sample}.cram",
    log:
        "results/processed_alignment/cram/{sample}.cram.log",
    threads: 2
    params:
        extra=lambda wildcards, input: f"-C -T {input.fa}",  # optional params string
        region="",  # optional region string
    wrapper:
        "v8.1.1/bio/samtools/view"


rule index_cram:
    input:
        rules.bam_to_cram.output,
    output:
        "results/processed_alignment/cram/{sample}.cram.crai",
    log:
        "results/processed_alignment/cram/{sample}_index.log",
    threads: 4  # This value - 1 will be sent to -@
    params:
        extra="",  # optional params string
    wrapper:
        "v8.1.1/bio/samtools/index"
