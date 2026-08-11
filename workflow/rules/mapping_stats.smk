rule gffread_gff:
    input:
        fasta=rules.get_genome.output.fasta,
        annotation=rules.get_genome.output.gff,
    output:
        records="results/get_genome/genome.bed",
    log:
        "results/get_genome/gffread.log",
    threads: 1
    params:
        extra=config["mapping_stats"]["gffread"]["extra"],
    message:
        "convert genome annotation from GFF to BED format"
    wrapper:
        "v5.0.0/bio/gffread"


rule rseqc_infer_experiment:
    input:
        aln=get_processed_bam,
        refgene="results/get_genome/genome.bed",
    output:
        "results/rseqc/infer_experiment/{sample}.txt",
    log:
        "results/rseqc/infer_experiment/{sample}.log",
    params:
        extra="--sample-size 10000",
    message:
        "infer experiment type from mapping to features"
    wrapper:
        "v4.7.5/bio/rseqc/infer_experiment"


rule rseqc_bam_stat:
    input:
        get_processed_bam,
    output:
        "results/rseqc/bam_stat/{sample}.txt",
    log:
        "results/rseqc/bam_stat/{sample}.log",
    threads: 2
    params:
        extra="--mapq 5",
    message:
        "collect mapping statistics using RSeQC"
    wrapper:
        "v5.0.0/bio/rseqc/bam_stat"


rule deeptools_coverage:
    input:
        bam=get_processed_bam,
        bai=get_processed_bam_index,
    output:
        "results/deeptools/coverage/{sample}.bw",
    log:
        "results/deeptools/coverage/{sample}.log",
    threads: 4
    params:
        effective_genome_size=config["mapping_stats"]["deeptools_coverage"][
            "genome_size"
        ],
        extra=config["mapping_stats"]["deeptools_coverage"]["extra"],
    message:
        "generate normalized coverage files using deeptools"
    wrapper:
        "v5.6.0/bio/deeptools/bamcoverage"
