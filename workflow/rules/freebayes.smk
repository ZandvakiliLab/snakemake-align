rule freebayes:
    input:
        alns=get_processed_bam,
        idxs=get_processed_bam_index,
        ref=rules.get_genome.output.fasta,
        refidx=rules.get_genome.output.fai,
    output:
        bcf="results/freebayes/call/{sample}.bcf",
    log:
        "results/freebayes/call/{sample}.log",
    threads: 4
    params:
        extra=config["variant_calling"]["freebayes"]["extra"],
    message:
        "call variants using freebayes"
    wrapper:
        "v9.4.1/bio/freebayes"
