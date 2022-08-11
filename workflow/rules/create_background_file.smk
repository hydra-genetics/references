# vim: syntax=python tabstop=4 expandtab
# coding: utf-8


__author__ = "Jonas Almlöf, Patrik Smeds"
__copyright__ = "Copyright 2021, Patrik Smeds, Jonas Almlöf"
__email__ = "jonas.almlöf@scilifelab.uu.se, patrik.smeds@scilifelab.uu.se"
__license__ = "GPL3"


rule create_background_file:
    input:
        gvcfs=get_gvcfs(units),
    output:
        background_file=temp("references/create_background_file/background_panel.tsv"),
    params:
        min_dp=config.get("create_artifact_file", {}).get("min_dp", 500),
        max_af=config.get("create_artifact_file", {}).get("max_af", 0.015),
    log:
        "references/create_background_file/create_background_file.log",
    conda:
        "../envs/create_background_file.yaml"
    container:
        config.get("create_background_file", {}).get("container", config["default_container"])
    script:
        "../scripts/create_background_file.py"
