rule gffread_gff:
    input:
        fasta="results/genome/genome.fasta",
        annotation="results/genome/genome.gff",
    output:
        records="results/genome/genome.bed",
    log:
        "results/genome/gffread.log",
    threads: 1
    params:
        extra=config["mapping_stats"]["gffread"]["extra"],
    message:
        "convert genome annotation from GFF to BED format"
    wrapper:
        "v5.0.0/bio/gffread"


rule rseqc_infer_experiment:
    input:
        aln=rules.filter_bam.output,
        refgene="results/genome/genome.bed",
    output:
        "results/rseqc/infer_experiment/{sample}_{genome}.txt",
    log:
        "results/rseqc/infer_experiment/{sample}_{genome}.log",
    params:
        extra="--sample-size 10000",
    message:
        "infer experiment type from mapping to features"
    wrapper:
        "v4.7.5/bio/rseqc/infer_experiment"


rule rseqc_bam_stat:
    input:
        rules.filter_bam.output,
    output:
        "results/rseqc/bam_stat/{sample}_{genome}.txt",
    log:
        "results/rseqc/bam_stat/{sample}_{genome}.log",
    threads: 2
    params:
        extra="--mapq 5",
    message:
        "collect mapping statistics using RSeQC"
    wrapper:
        "v5.0.0/bio/rseqc/bam_stat"

rule deeptools_coverage_without_normalization:
    # Good for sparse or empty BAM files
    input:
        bam=rules.filter_bam.output,
        bai=rules.samtools_index_processed.output,
    output:
        temp("results/deeptools/coverage/{sample}.{genome}.{strand}.nonorm.bw"),
    log:
        "results/deeptools/coverage/{sample}.{genome}.{strand}.nonorm.log",
    wildcard_constraints:
        strand="plus|minus",
        sample="|".join(samples.index),
    threads: 4
    params:
        effective_genome_size=lambda wc: config["mapping_stats"]["deeptools_coverage"][
            "genome_size"
        ][wc.genome],
        extra=lambda wc: (
            config["mapping_stats"]["deeptools_coverage"]["extra"]
            + " --filterRNAstrand {strand}".format(
                strand="forward" if wc.strand == "plus" else "reverse"
            )
        ),
    message:
        "generate {wildcards.strand}-strand coverage using deeptools for sample {wildcards.sample}"
    wrapper:
        "v5.6.0/bio/deeptools/bamcoverage"

rule deeptools_coverage:
    input:
        bam=rules.filter_bam.output,
        bai=rules.samtools_index_processed.output,
    output:
        temp("results/deeptools/coverage/{sample}.{genome}.{strand}.norm.bw"),
    log:
        "results/deeptools/coverage/{sample}.{genome}.{strand}.norm.log",
    wildcard_constraints:
        strand="plus|minus",
        sample="|".join(samples.index),
    threads: 4
    params:
        effective_genome_size=lambda wc: config["mapping_stats"]["deeptools_coverage"][
            "genome_size"
        ][wc.genome],
        extra=lambda wc: (
            config["mapping_stats"]["deeptools_coverage"]["extra"]
            + "--normalizeUsing {norm}".format(
                norm=config["mapping_stats"]["deeptools_coverage"]["normalizeUsing"]
            )
            + " --filterRNAstrand {strand}".format(
                strand="forward" if wc.strand == "plus" else "reverse"
            )
        ),
    message:
        "generate normalized {wildcards.strand}-strand coverage using deeptools for sample {wildcards.sample}"
    wrapper:
        "v5.6.0/bio/deeptools/bamcoverage"

checkpoint samtools_count:
    input:
        rules.filter_bam.output,
    output:
        "results/processed_alignment/filtered_bam/{genome}/{sample}.count",
    log:
        "results/processed_alignment/filtered_bam/{genome}/{sample}.count.log",
    wildcard_constraints:
        strand="plus|minus",
        sample="|".join(samples.index),
    params:
        extra="-c",
    message:
        "count number of reads aligning to {wildcards.genome} in sample {wildcards.sample}"
    conda:
        "../envs/samtools.yml"
    shell:
        "samtools view -c {input} > {output} 2> {log}"

rule bigwig:
    input:
        choose_bamCoverage_method
    output:
        "results/deeptools/coverage/{sample}.{genome}.{strand}.bw"
    log:
        "results/deeptools/coverage/{sample}.{genome}.{strand}.move.log"
    conda:
        "../envs/basic.yml"
    shell:
        """
        mv {input} {output}
        """