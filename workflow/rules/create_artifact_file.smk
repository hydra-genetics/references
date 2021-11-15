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
        bed=temp("references/cnvkit_create_targets/cnvkit_manifest.target.bed"),
    log:
        "references/cnvkit_create_targets/cnvkit_create_targets.log",
    conda:
        "../envs/cnvkit_panel_of_normal.yaml"
    container:
        config.get("cnvkit_create_targets", {}).get("container", config["default_container"])
    shell:
        "(cnvkit.py target --split {input.bed} -o {output.bed}) &> {log}"
