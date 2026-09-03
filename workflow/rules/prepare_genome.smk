# ─────────────────────────────────────────────────────────────────────────────
# 1. Get Genome
# ─────────────────────────────────────────────────────────────────────────────


rule efetch_accession:
    output:
        "results/ncbi/{ref}_accession.txt",
    log:
        "results/ncbi/{ref}_accession.log",
    conda:
        "../envs/efetch.yml"
    params:
        genbank_id=lambda wc: REFS[wc.ref]["datasets"]["genbank"],
    message:
        "get genome assembly accession from genbank id"
    shell:
        """
        esearch -db nuccore -query {params.genbank_id} 2>{log} \
            | elink -target assembly \
            | esummary | xtract -pattern DocumentSummary -element AssemblyAccession \
            >{output}

        # Verify we actually got something
        if [ ! -s {output} ]; then
            echo "ERROR: Could not resolve assembly accession for {params.genbank_id}" >&2
            exit 1
        fi
        """


rule ncbi_datasets:
    input:
        branch(
            lambda wc: REFS[wc.ref]["datasets"]["assembly"],
            then=[],
            otherwise=rules.efetch_accession.output,
        ),
    output:
        fasta="results/ncbi/{ref}.fasta",
        gff="results/ncbi/{ref}.gff3",
    log:
        "results/ncbi/{ref}.log",
    conda:
        "../envs/ncbi_datasets.yml"
    params:
        accession=lambda wc: REFS[wc.ref]["datasets"]["assembly"],
        extra=lambda wc: REFS[wc.ref]["datasets"]["extra"] or "",
    message:
        "download genome genome from NCBI"
    shell:
        """
        # Use a local scratch directory
        WORKDIR=$(mktemp -d)
        trap "rm -rf $WORKDIR" EXIT

        # Determine genome assembly accession
        if [[ -n "{params.accession}" ]]; then
            ACCESSION="{params.accession}"
        else
            ACCESSION=$(tail -1 {input})
        fi

        echo "Using accession: $ACCESSION" >{log}

        # Download genome
        datasets download genome accession "$ACCESSION" \
            --include genome,gff3 {params.extra} \
            --filename "$WORKDIR/genome.zip" >>{log} 2>&1

        # Unzip
        unzip -o "$WORKDIR/genome.zip" -d "$WORKDIR/genome" >>{log} 2>&1

        # Collect FASTA
        find "$WORKDIR/genome" -name "*.fna" | sort | xargs cat >{output.fasta}

        # Collect GFF3
        find "$WORKDIR/genome" -name "*.gff" | sort | xargs cat >{output.gff}
        """


rule get_genome:
    input:
        branch(
            lambda wc: REFS[wc.ref]["source"] == "local",
            then=lambda wc: REFS[wc.ref]["local"]["fasta"],
            otherwise=rules.ncbi_datasets.output.fasta,
        ),
    output:
        "results/genome/{ref}.fasta",
    log:
        "results/genome/{ref}_fasta.log",
    wildcard_constraints:
        ref="|".join(REF_KEYS),
    conda:
        "../envs/basic.yml"
    message:
        "obtaining fasta files"
    shell:
        """
        input_file=$(realpath {input})
        ln -s $input_file {output}
        echo 'made symbolic link from {input} to {output}' >{log}
        """


rule get_gff:
    input:
        branch(
            lambda wc: REFS[wc.ref]["source"] == "local",
            then=lambda wc: REFS[wc.ref]["local"]["gff"],
            otherwise=rules.ncbi_datasets.output.gff,
        ),
    output:
        "results/genome/{ref}.gff",
    log:
        "results/genome/{ref}_gff.log",
    wildcard_constraints:
        ref="|".join(REF_KEYS),
    conda:
        "../envs/basic.yml"
    shell:
        """
        input_file=$(realpath -e {input})
        ln -sf $input_file {output}
        echo 'made symbolic link from {input} to {output}' >{log}
        """


# ─────────────────────────────────────────────────────────────────────────────
# 2. Add contig-name prefixes to each individual FASTA & GFFs
# ─────────────────────────────────────────────────────────────────────────────


rule add_prefix_to_fasta:
    """Prepend '<prefix>_' to every contig name in a reference FASTA."""
    input:
        rules.get_genome.output,
    output:
        "results/genome/{ref}_prefixed.fasta",
    log:
        "results/genome/{ref}_prefixed.fasta.log",
    conda:
        "../envs/basic.yml"
    params:
        prefix=lambda wc: REFS[wc.ref]["prefix"],
    shell:
        """
        awk -v pfx="{params.prefix}" '
            /^>/ {{ print ">" pfx "_" substr($0, 2) }}
            !/^>/ {{ print }}
        ' {input} >{output} 2>{log}
        """


rule add_prefix_to_gff:
    """
    Add a prefix to every contig name (column 1) in a GFF file.
    Handles standard GFF2/GFF3 format; skips comment lines.
    """
    input:
        rules.get_gff.output,
    output:
        "results/genome/{ref}_prefixed.gff",
    log:
        "results/genome/{ref}_prefixed.gff.log",
    conda:
        "../envs/basic.yml"
    params:
        prefix=lambda wc: REFS[wc.ref]["prefix"],
    shell:
        """
        awk -v pfx="{params.prefix}" '
            /^#/ {{ print; next }}
            NF   {{
                idx = index($0, "\\t")
                print pfx "_" substr($0, 1, idx-1) substr($0, idx)
            }}
        ' {input} >{output} 2>{log}
        """


# ─────────────────────────────────────────────────────────────────────────────
# 3. Concatenate Genomes
# ─────────────────────────────────────────────────────────────────────────────


rule concatenate_reference:
    """Concatenate an arbitrary number of prefixed reference FASTAs."""
    input:
        get_genome_files,
    output:
        "results/genome/genome.{filetype}",
    log:
        "results/genome/genome.{filetype}.log",
    conda:
        "../envs/basic.yml"
    message:
        "concatenate fasta files"
    shell:
        "cat {input} > {output} 2> {log}"
