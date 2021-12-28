# vim: syntax=python tabstop=4 expandtab
# coding: utf-8


__author__ = "Jonas Almlöf, Patrik Smeds"
__copyright__ = "Copyright 2021, Patrik Smeds, Jonas Almlöf"
__email__ = "jonas.almlöf@scilifelab.uu.se, patrik.smeds@scilifelab.uu.se"
__license__ = "GPL3"

import os


rule svdb_build:
    input:
        cnv_vcfs=get_cnv_vcfs(units),
    output:
        svdb_db="references/svdb_build/svdb_cnv.db",
    params:
        prefix=lambda wildcards, output: os.path.splitext(output[0])[0],
    log:
        "references/svdb_build/svdb_build.log",
    conda:
        "../envs/svdb_cnv_db.yaml"
    container:
        config.get("svdb_build", {}).get("container", config["default_container"])
    shell:
        "svdb --build --prefix {params.prefix} --files {input.cnv_vcfs}"


rule svdb_export:
    input:
        svdb_db="references/svdb_build/svdb_cnv.db",
    output:
        svdb_vcf="references/svdb_export/svdb_cnv.vcf",
    params:
        prefix=lambda wildcards, output: os.path.splitext(output[0])[0],
        extra=config.get("svdb_export", {}).get("extra", ""),
    log:
        "references/svdb_export/svdb_export.log",
    conda:
        "../envs/svdb_cnv_db.yaml"
    container:
        config.get("svdb_export", {}).get("container", config["default_container"])
    shell:
        "svdb --export --db {input.svdb_db} --prefix {params.prefix} {params.extra}"
