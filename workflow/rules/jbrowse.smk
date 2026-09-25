# ─────────────────────────────────────────────────────────────────────────────
# 1. Add Assembly
# ─────────────────────────────────────────────────────────────────────────────


rule faToTwoBit_fa:
    input:
        "results/genome/{genome}.fasta",
    output:
        "results/genome/{genome}.2bit",
    log:
        "results/genome/{genome}.fa_to_2bit.log",
    wrapper:
        "v7.1.0/bio/ucsc/faToTwoBit"


rule jbrowse_add_assembly:
    input:
        fa="results/genome/{genome}.2bit",
    output:
        config=temp("jbrowse/config_{genome}_assembly.json"),
    log:
        "results/jbrowse/add_assembly_{genome}.log",
    conda:
        "../envs/jbrowse.yml"
    resources:
        file_lock=1,
    params:
        clean_path=lambda wc, input: strip_prefix(input.fa),
        url_prefix=config["jbrowse"]["path_prefix"],
        extra=lambda wc: config["jbrowse"]["add_assembly"][wc.genome],
    message:
        "add {wildcards.genome} assemblies to jbrowse"
    shell:
        """
        (
            # Add to jbrowse
            path_to_asm="{params.url_prefix}{params.clean_path}"
            jbrowse add-assembly "$path_to_asm" \
                --type twoBit \
                --target {output.config} \
                ${params.extra}
        ) >{log} 2>&1
        """


# ─────────────────────────────────────────────────────────────────────────────
# 2. Add Annotation
# ─────────────────────────────────────────────────────────────────────────────


rule sort_gff:
    input:
        gff="results/genome/{genome}.gff",
    output:
        gff="results/genome/{genome}.sorted.gff.gz",
    log:
        "results/genome/{genome}_sort_gff.log",
    conda:
        "../envs/jbrowse.yml"
    message:
        "sort gff3"
    shell:
        """
        jbrowse sort-gff {input.gff} | bgzip >{output.gff}
        """


rule index_gff:
    input:
        "results/genome/{genome}.sorted.gff.gz",
    output:
        "results/genome/{genome}.sorted.gff.gz.tbi",
    log:
        "results/genome/{genome}_index_gff.log",
    conda:
        "../envs/jbrowse.yml"
    message:
        "index gff3"
    shell:
        """
        tabix {input}
        """


rule jbrowse_add_anno:
    input:
        gff="results/genome/{genome}.sorted.gff.gz",
        tbi="results/genome/{genome}.sorted.gff.gz.tbi",
        config="jbrowse/config_{genome}_assembly.json",
    output:
        config=temp("jbrowse/config_{genome}_anno.json"),
    log:
        "results/jbrowse/add_anno_{genome}.log",
    conda:
        "../envs/jbrowse.yml"
    resources:
        file_lock=1,
    params:
        clean_path_gff=lambda wc, input: strip_prefix(input.gff),
        clean_path_tbi=lambda wc, input: strip_prefix(input.tbi),
        url_prefix=config["jbrowse"]["path_prefix"],
        extra=lambda wc: config["jbrowse"]["add_anno"][wc.genome],
    message:
        "add {wildcards.genome} annotations to jbrowse"
    shell:
        """
        (

            # Get url/path to annotation files
            path_to_gff="{params.url_prefix}{params.clean_path_gff}"
            path_to_tbi="{params.url_prefix}{params.clean_path_tbi}"

            # add to jbrowse
            cp {input.config} {output.config}
            jbrowse add-track "$path_to_gff" \
                --indexFile "$path_to_tbi" \
                --target {output.config} \
                --assemblyNames {wildcards.genome} \
                {params.extra}
        ) >{log} 2>&1
        """


# ─────────────────────────────────────────────────────────────────────────────
# 3. Add BigWigs
# ─────────────────────────────────────────────────────────────────────────────


rule jbrowse_add_bw:
    input:
        config="jbrowse/config_{genome}_anno.json",
        plus_bw=expand(
            "results/deeptools/coverage/{sample}.{{genome}}.plus.bw",
            sample=samples.index,
        ),
        minus_bw=expand(
            "results/deeptools/coverage/{sample}.{{genome}}.minus.bw",
            sample=samples.index,
        ),
    output:
        config=temp("jbrowse/config_{genome}_bw.json"),
    log:
        "results/jbrowse/add_bw_{genome}.log",
    conda:
        "../envs/jbrowse.yml"
    resources:
        file_lock=1,
    params:
        clean_path_plus=lambda wc, input: strip_prefix(input.plus_bw),
        clean_path_minus=lambda wc, input: strip_prefix(input.minus_bw),
        url_prefix=config["jbrowse"]["path_prefix"],
        extra=config["jbrowse"]["add_bw"]["extra"],
    message:
        "add BigWigs to {wildcards.genome} config.json"
    shell:
        """
        (
            cp {input.config} {output.config}

            for i in {params.clean_path_plus}; do

                # Get url/path to annotation files
                path_to_bw="{params.url_prefix}$i"

                # Add to jbrowse
                jbrowse add-track $path_to_bw \
                    --target {output.config} \
                    --name "${{i##*/}}" \
                    --assemblyNames {wildcards.genome} \
                    {params.extra}
            done

            for i in {params.clean_path_minus}; do

                # Get url/path to annotation files
                path_to_bw="{params.url_prefix}$i"

                # Add to jbrowse
                jbrowse add-track $path_to_bw \
                    --target {output.config} \
                    --name "${{i##*/}}" \
                    --assemblyNames {wildcards.genome} \
                    --config '{{"displays":[{{"type":"LinearWiggleDisplay","displayId":"my_bw-LinearWiggleDisplay","inverted":true}}]}}' \
                    {params.extra}
            done

        ) >{log} 2>&1
        """


# ─────────────────────────────────────────────────────────────────────────────
# 4. Add cram files
# ─────────────────────────────────────────────────────────────────────────────


rule jbrowse_add_cram:
    input:
        config="jbrowse/config_{genome}_bw.json",
        cram=expand(
            "results/processed_alignment/cram/{sample}_{{genome}}.cram",
            sample=samples.index,
        ),
    output:
        config=temp("jbrowse/config_{genome}_cram.json"),
    log:
        "results/jbrowse/add_cram_{genome}.log",
    conda:
        "../envs/jbrowse.yml"
    resources:
        file_lock=1,
    params:
        clean_path=lambda wc, input: strip_prefix(input.cram),
        url_prefix=config["jbrowse"]["path_prefix"],
        extra=config["jbrowse"]["add_cram"]["extra"],
    message:
        "add plus cram tracks to {wildcards.genome} config.json"
    shell:
        """
        (
            cp {input.config} {output.config}
            for i in {input.cram}; do

                # Get url/path to annotation files
                path_to_cram="{params.url_prefix}$i"

                # Add to jbrowse
                jbrowse add-track $path_to_cram \
                    --indexFile $i.crai \
                    --target {output.config} \
                    --name "${{i##*/}}" \
                    --assemblyNames {wildcards.genome} \
                    --config '{{"displays":[{{"type":"LinearPileupDisplay", "colorBySetting": {{"type": "strand"}}}}]}}' \
                    {params.extra}
            done
        ) >{log} 2>&1
        """


# ─────────────────────────────────────────────────────────────────────────────
# 5. Merge config.jsons
# ─────────────────────────────────────────────────────────────────────────────


rule jbrowse_merge_jsons:
    input:
        config_files=expand(
            "jbrowse/config_{genome}_cram.json",
            genome=list(config["jbrowse"]["add_assembly"].keys()),
        ),
    output:
        "jbrowse/config.json",
    log:
        "results/jbrowse/merge_config_json.log",
    conda:
        "../envs/jbrowse.yml"
    message:
        "merge config files"
    shell:
        """
        jq -s '
            {{
                assemblies:                  map(.assemblies)                  | add,
                tracks:                      map(.tracks)                      | add,
                connections:                 map(.connections // [])           | add,
                aggregateTextSearchAdapters: map(.aggregateTextSearchAdapters // []) | add,
                configuration:               map(.configuration)               | add,
                defaultSession:              .[0].defaultSession
            }}
        ' {input.config_files} >{output}
        """
