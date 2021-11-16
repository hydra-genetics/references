# vim: syntax=python tabstop=4 expandtab
# coding: utf-8


__author__ = "Jonas Almlöf, Patrik Smeds"
__copyright__ = "Copyright 2021, Patrik Smeds, Jonas Almlöf"
__email__ = "jonas.almlöf@scilifelab.uu.se, patrik.smeds@scilifelab.uu.se"
__license__ = "GPL3"


rule create_artifact_file:
    input:
        vcfs=get_vcfs(units),
    output:
        artifact_panel="references/create_artifact_file/artifact_panel.tsv",
    params:
        callers=["vardict", "mutect2", "freebayes", "varscan"],
    log:
        "references/create_artifact_file/create_artifact_file.log",
    conda:
        "../envs/create_artifact_file.yaml"
    container:
        config.get("cnvkit_create_targets", {}).get("container", config["default_container"])
    script:
        "../scripts/create_artifact_file.py"
