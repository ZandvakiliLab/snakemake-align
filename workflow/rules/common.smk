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
# GENOME-RELATED FUNCTIONS
###############################

REFS = config["genome"]  # dict of all reference blocks
REF_KEYS = list(REFS.keys())  # e.g. ["host", "virus", "spike"]
SPIKE_KEY = "spike"  # which reference is the spike-in


def get_genome_files(wildcards):
    files = []
    for genome, vals in REFS.items():
        if vals["prefix"]:
            files.append(f"results/genome/{genome}_prefix.{wildcards.filetype}")
        else:
            files.append(f"results/genome/{genome}.{wildcards.filetype}")
    return files

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


# get processed fastq files (after fastp or umi_tools)
def get_processed_fastq(wildcards, regex=None):
    processed_fastq = expand(
        "{dir}/{{sample}}_{read}.fastq.gz",
        read=["read1", "read2"] if is_paired_end() else ["read1"],
        dir=config["mapping"]["fastq_dir"],
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


def get_bam_2(wildcards):
    if (
        config["mapping_postprocessing"]["deduplication"]["enabled"]
        and config["mapping_postprocessing"]["deduplication"]["tool"] == "samtools"
    ):
        return rules.samtools_fixmate.output
    else:
        return get_bam(wildcards)


def get_processed_bam(wildcards):
    if config["mapping_postprocessing"]["deduplication"]["enabled"]:
        return expand(
            "results/processed_alignment/dedup/{tool}/{sample}.bam",
            sample=wildcards.sample,
            tool=config["mapping_postprocessing"]["deduplication"]["tool"],
        )
    else:
        return rules.samtools_sort.output


def get_processed_bam_index(wildcards):
    if config["mapping_postprocessing"]["deduplication"]["enabled"]:
        return expand(
            "results/processed_alignment/dedup/{tool}/{sample}.bam.bai",
            sample=wildcards.sample,
            tool=config["mapping_postprocessing"]["deduplication"]["tool"],
        )
    else:
        return rules.samtools_index.output


####################
# MULTIQC FUNCTION
####################


# get input for multiqc
def get_multiqc_input(wildcards):
    result = []
    result += expand(
        "results/{tool}/align/{sample}/mapped.bam",
        sample=samples.index,
        tool=config["mapping"]["tool"],
    )
    if config["mapping_postprocessing"]["deduplication"]["enabled"]:
        if config["mapping_postprocessing"]["deduplication"]["tool"] == "samtools":
            result += expand(
                "results/processed_alignment/dedup/samtools/{sample}_markdup.json",
                sample=samples.index,
            )
        else:
            result += expand(
                "results/processed_alignment/dedup/{tool}/{sample}.log",
                sample=samples.index,
                tool=config["mapping_postprocessing"]["deduplication"]["tool"],
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
    return result
