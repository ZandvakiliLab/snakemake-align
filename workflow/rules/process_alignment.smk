rule samtools_sort:
    input:
        get_bam,
    output:
        "results/processed_alignment/sort/{sample}.bam",
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
        "results/processed_alignment/sort/{sample}.bam.bai",
    log:
        "results/processed_alignment/sort/{sample}_index.log",
    threads: 2
    params:
        extra=config["mapping_postprocessing"]["samtools_index"]["extra"],
    message:
        "index reads"
    wrapper:
        "v9.4.1/bio/samtools/index"


rule umi_tools_dedup:
    input:
        bam=rules.samtools_sort.output,
        bai=rules.samtools_index.output,
    output:
        "results/processed_alignment/dedup/{sample}.bam",
    log:
        "results/processed_alignment/dedup/{sample}.log",
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


rule samtools_index_dedup:
    input:
        rules.umi_tools_dedup.output,
    output:
        "results/processed_alignment/dedup/{sample}.bam.bai",
    log:
        "results/processed_alignment/dedup/{sample}_index.log",
    threads: 2
    params:
        extra=config["mapping_postprocessing"]["samtools_index"]["extra"],
    message:
        "index reads"
    wrapper:
        "v9.4.1/bio/samtools/index"
