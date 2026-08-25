rule multiqc:
    input:
        get_multiqc_input,
    output:
        report="results/multiqc/multiqc_report.html",
    log:
        "results/multiqc/multiqc.log",
    params:
        extra=config["qc"]["multiqc"]["extra"],
        use_input_files_only=True,
    message:
        "generating MultiQC report for seq data"
    wrapper:
        "v8.1.1/bio/multiqc"
