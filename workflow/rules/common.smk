# import basic packages
import pandas as pd
from snakemake.utils import validate
from pathlib import Path

# read sample sheet
samples = (
    pd.read_csv(config["samplesheet"], sep="\t", dtype={"sample": str})
    .set_index("sample", drop=False)
    .sort_index()
)


# validate sample sheet and config file
validate(samples, schema="../../config/schemas/samples.schema.yml")
validate(config, schema="../../config/schemas/config.schema.yml")

###############################
# FASTQ-RELATED FUNCTIONS
###############################


# determine input type
def is_paired_end():
    if samples["read2"].isna().all():
        return False
    elif samples["read2"].notna().all():
        return True
    else:
        raise ValueError(
            f"Some samples seem to have a read2 fastq file, while others have only a "
            + "read1 fastq file. \nYou may not mix single-end and paired-end samples."
        )


# get fastq files
def get_fastq(wildcards):
    file = Path(samples.loc[wildcards["sample"]][wildcards["read"]])
    if file.is_absolute():
        return file
    else:
        input_dir = Path.absolute(Path.cwd())
        return input_dir / file


# get pairs of fastq files for fastp
def get_fastq_pairs(wildcards):
    return expand(
        "results/{folder}/{sample}_{read}.fastq.gz",
        folder=(
            "umi_tools/extract"
            if config["processing"]["umi_extraction"]["enabled"]
            else "get_fastq"
        ),
        sample=wildcards.sample,
        read=["read1", "read2"] if is_paired_end() else ["read1"],
    )


# get processed fastq files (after fastp or umi_tools)
def get_processed_fastq(wildcards, regex=None):
    processed_fastq = expand(
        "results/{tool}/{{sample}}_{read}.fastq.gz",
        read=["read1", "read2"] if is_paired_end() else ["read1"],
        tool=config["processing"]["tool"],
    )

    if regex is None:
        return processed_fastq
    else:
        return [s for s in processed_fastq if re.search(regex, s)]


###############################
# ALIGNMENT-RELATED FUNCTIONS
###############################


# get bam files
def get_bam(wildcards):
    return expand(
        "results/{tool}/align/{sample}/mapped.bam",
        sample=wildcards.sample,
        tool=config["mapping"]["tool"],
    )


###############################
# ALIGNMENT PROCESSING-RELATED FUNCTIONS
###############################


def get_processed_bam(wildcards):
    if config["mapping_postprocessing"]["deduplication"]["enabled"]:
        return rules.umi_tools_dedup.output
    else:
        return rules.samtools_sort.output


def get_processed_bam_index(wildcards):
    if config["mapping_postprocessing"]["deduplication"]["enabled"]:
        return rules.samtools_index_dedup.output
    else:
        return rules.samtools_index.output


###############################
# VARIANT-RELATED FUNCTIONS
###############################


# get variant files to make consensus
def get_variants(wildcards):
    return expand(
        "results/{caller}/effect/{sample}{ext}.gz",
        caller=wildcards.caller,
        sample=samples.index,
        ext=(
            "_vep.vcf"
            if config["variant_annotation"]["tool"] == "vep"
            else "_snpeff.vcf"
        ),
    )


####################
# MULTIQC FUNCTION
####################


# get input for multiqc
def get_multiqc_input(wildcards):
    result = []
    result += expand(
        "results/fastqc/{sample}_{read}_fastqc.{ext}",
        sample=samples.index,
        read=["read1", "read2"] if is_paired_end() else ["read1"],
        ext=["html", "zip"],
    )
    result += expand(
        "results/{tool}/align/{sample}/mapped.bam",
        sample=samples.index,
        tool=config["mapping"]["tool"],
    )
    result += expand(
        "results/fastp/{sample}.json",
        sample=samples.index,
    )
    result += expand(
        "results/rseqc/{tool}/{sample}.txt",
        sample=samples.index,
        tool=["infer_experiment", "bam_stat"],
    )
    result += expand(
        "results/deeptools/coverage/{sample}.bw",
        sample=samples.index,
    )
    result += expand(
        "results/{caller}/call/{sample}{ext}",
        sample=samples.index,
        caller=config["variant_calling"]["tool"],
        ext=["_stats.txt", "_all.vcf", "_variants.vcf"],
    )
    result += expand(
        "results/{caller}/effect/{sample}{ext}",
        sample=samples.index,
        caller=config["variant_calling"]["tool"],
        ext=(
            ["_vep.vcf", "_vep.html"]
            if config["variant_annotation"]["tool"] == "vep"
            else ["_snpeff.vcf", "_snpeff.csv"]
        ),
    )
    return result
