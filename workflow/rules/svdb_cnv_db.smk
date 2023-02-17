__author__ = "Jonas Almlöf, Patrik Smeds"
__copyright__ = "Copyright 2021, Patrik Smeds, Jonas Almlöf"
__email__ = "jonas.almlöf@scilifelab.uu.se, patrik.smeds@scilifelab.uu.se"
__license__ = "GPL3"

import os


rule svdb_build:
    input:
        cnv_vcfs=get_cnv_vcfs(units),
    output:
        svdb_db=temp("references/svdb_build/svdb_cnv.db"),
    params:
        prefix=lambda wildcards, output: os.path.splitext(output[0])[0],
    log:
        "references/svdb_build/svdb_cnv.db.log",
    benchmark:
        repeat(
            "references/svdb_build/svdb_cnv.db.benchmark.tsv",
            config.get("svdb_build", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("svdb_build", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("svdb_build", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("svdb_build", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("svdb_build", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("svdb_build", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("svdb_build", {}).get("time", config["default_resources"]["time"]),
    conda:
        "../envs/svdb_cnv_db.yaml"
    container:
        config.get("svdb_build", {}).get("container", config["default_container"])
    message:
        "{rule}: build svdb DB from vcfs"
    shell:
        "svdb --build --prefix {params.prefix} --files {input.cnv_vcfs}"


rule svdb_export:
    input:
        svdb_db="references/svdb_build/svdb_cnv.db",
    output:
        svdb_vcf=temp("references/svdb_export/svdb_cnv.vcf"),
    params:
        prefix=lambda wildcards, output: os.path.splitext(output[0])[0],
        extra=config.get("svdb_export", {}).get("extra", ""),
    log:
        "references/svdb_export/svdb_cnv.vcf.log",
    benchmark:
        repeat(
            "references/svdb_export/svdb_cnv.vcf.benchmark.tsv",
            config.get("svdb_export", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("svdb_export", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("svdb_export", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("svdb_export", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("svdb_export", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("svdb_export", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("svdb_export", {}).get("time", config["default_resources"]["time"]),
    conda:
        "../envs/svdb_cnv_db.yaml"
    container:
        config.get("svdb_export", {}).get("container", config["default_container"])
    message:
        "{rule}: export svdb db to vcf"
    shell:
        "svdb --export --db {input.svdb_db} --prefix {params.prefix} {params.extra}"
