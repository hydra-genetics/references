# vim: syntax=python tabstop=4 expandtab
# coding: utf-8


__author__ = "Jonas Almlöf, Patrik Smeds"
__copyright__ = "Copyright 2021, Patrik Smeds, Jonas Almlöf"
__email__ = "jonas.almlöf@scilifelab.uu.se, patrik.smeds@scilifelab.uu.se"
__license__ = "GPL3"


rule cnvkit_create_targets:
    input:
        bed=config["reference"]["design_bedfile"],
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


rule cnvkit_create_anti_targets:
    input:
        bed="references/cnvkit_create_targets/cnvkit_manifest.target.bed",
    output:
        bed=temp("references/cnvkit_create_anti_targets/cnvkit_manifest.antitarget.bed"),
    log:
        "references/cnvkit_create_anti_targets/cnvkit_create_anti_targets.log",
    conda:
        "../envs/cnvkit_panel_of_normal.yaml"
    container:
        config.get("cnvkit_create_anti_targets", {}).get("container", config["default_container"])
    shell:
        "(cnvkit.py antitarget {input.bed} -o {output.bed}) &> {log}"


rule cnvkit_build_normal_reference:
    input:
        bams=get_bams(units),
        target="references/cnvkit_create_targets/cnvkit_manifest.target.bed",
        antitarget="references/cnvkit_create_anti_targets/cnvkit_manifest.antitarget.bed",
        ref=config["reference"]["fasta"],
        mappability=config["reference"]["mappability"],
    output:
        PoN="references/cnvkit_build_normal_reference/cnvkit.PoN.cnn",
        tmp_bed=temp("cnvkit_manifest.target.target.bed"),
        tmp_target_cov=temp(
            ["%s_%s.targetcoverage.cnn" % (sample, t) for sample in get_samples(samples) for t in get_unit_types(units, sample)]
        ),
        tmp_antitarget_cov=temp(
            [
                "%s_%s.antitargetcoverage.cnn" % (sample, t)
                for sample in get_samples(samples)
                for t in get_unit_types(units, sample)
            ]
        ),
    params:
        extra=config.get("cnvkit_build_normal_reference", {}).get("extra", ""),
    log:
        "references/cnvkit_build_normal_reference/cnvkit_build_normal_reference.log",
    threads: 4
    conda:
        "../envs/cnvkit_panel_of_normal.yaml"
    container:
        config.get("cnvkit_build_normal_reference", {}).get("container", config["default_container"])
    shell:
        "(cnvkit.py batch "
        " {params.extra} "
        "-n {input.bams} "
        "-m hybrid "
        "--output-reference {output.PoN} "
        "-t {input.target} "
        "-f {input.ref} "
        "-a {input.antitarget} "
        "-g {input.mappability} "
        "-p {threads}) &> {log}"
