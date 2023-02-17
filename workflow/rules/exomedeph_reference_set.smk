__author__ = "Padraic Corcoran"
__copyright__ = "Copyright 2023, Padraic Corcoran"
__email__ = "padraic.corcoran@scilifelab.uu.se"
__license__ = "GPL3"


rule exomedepth_bam_list:
    input:
        bam_list=get_bams(units),
    output:
        bam_list_file="references/exomedepth_bam_list/bam_files.list",
    log:
        "references/exomedepth_bam_list/bam_files.list.log",
    benchmark:
        repeat(
            "references/exomedepth_bam_list/bam_files.list.benchmark.tsv",
            config.get("exomedepth_bam_list", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("exomedepth_bam_list", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("exomedepth_bam_list", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("exomedepth_bam_list", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("exomedepth_bam_list", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("exomedepth_bam_list", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("exomedepth_bam_list", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("exomedepth_bam_list", {}).get("container", config["default_container"])
    conda:
        "../envs/exomedepth_reference.yaml"
    message:
        "{rule}: make a bam list file for purecn coverage"
    shell:
        "(for val in {input.bam_list}; do echo $val >> {output.bam_list_file}; done) &> {log}"


rule exomedepth_reference:
    input:
        bam_list_file="references/exomedepth_bam_list/bam_files.list",
    output:
        reference="references/exomedepth_reference/RefCount.mat"
    params:
        bed=config.get("exomedepth_reference", {}).get("exons_bed", ""),
        extra=config.get("exomedepth_reference", {}).get("extra", ""),
    log:
        "references/exomedepth_reference/bam_files.list.output.log",
    benchmark:
        repeat(
            "references/exomedepth_reference/bam_files.list.benchmark.tsv",
            config.get("exomedepth_reference", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("exomedepth_reference", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("exomedepth_reference", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("exomedepth_reference", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("exomedepth_reference", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("exomedepth_reference", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("exomedepth_reference", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("exomedepth_reference", {}).get("container", config["default_container"])
    conda:
        "../envs/exomedepth_reference.yaml"
    message:
        "{rule}: calculate read counts for samples in {input.bam_list_file}"
    script:
        "exomedepth_reference.R"
