rule minimap2_index:
    input:
        target="results/genome/genome.fasta",
    output:
        index="results/minimap2/index/genome.mmi",
    log:
        "results/minimap2/index/genome.log",
    threads: 1
    params:
        extra=config["mapping"]["minimap2"]["index"],
    wrapper:
        "v7.1.0/bio/minimap2/index"


rule minimap2_align:
    input:
        target=rules.minimap2_index.output.index,
        query=get_processed_fastq,
    output:
        "results/minimap2/align/{sample}/mapped.bam",
    log:
        "results/minimap2/align/{sample}/mapped.log",
    threads: 8
    params:
        extra=config["mapping"]["minimap2"]["extra"],
        sorting=config["mapping"]["minimap2"]["sorting"],
        sort_extra=config["mapping"]["minimap2"]["sort_extra"],
    wrapper:
        "v9.4.1/bio/minimap2/aligner"
