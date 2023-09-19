__author__ = "Jonas Almlöf, Patrik Smeds"
__copyright__ = "Copyright 2021, Patrik Smeds, Jonas Almlöf"
__email__ = "jonas.almlöf@scilifelab.uu.se, patrik.smeds@scilifelab.uu.se"
__license__ = "GPL3"


rule cnvkit_create_targets:
    input:
        bed=config.get("reference", {}).get("design_bed", ""),
    output:
        bed=temp("references/cnvkit_create_targets/cnvkit_manifest.target.bed"),
    params:
        extra=config.get("cnvkit_create_targets", {}).get("extra", ""),
    log:
        "references/cnvkit_create_targets/cnvkit_manifest.target.bed.log",
    benchmark:
        repeat(
            "references/cnvkit_create_targets/cnvkit_manifest.target.bed.benchmark.tsv",
            config.get("cnvkit_create_targets", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("cnvkit_create_targets", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("cnvkit_create_targets", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("cnvkit_create_targets", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("cnvkit_create_targets", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("cnvkit_create_targets", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("cnvkit_create_targets", {}).get("time", config["default_resources"]["time"]),
    conda:
        "../envs/cnvkit_panel_of_normal.yaml"
    container:
        config.get("cnvkit_create_targets", {}).get("container", config["default_container"])
    message:
        "{rule}: create target bed file for cnvkit from design bed"
    shell:
        "(cnvkit.py target --split {input.bed} -o {output.bed} {params.extra}) &> {log}"


rule cnvkit_create_anti_targets:
    input:
        bed="references/cnvkit_create_targets/cnvkit_manifest.target.bed",
    output:
        bed=temp("references/cnvkit_create_anti_targets/cnvkit_manifest.antitarget.bed"),
    params:
        extra=config.get("cnvkit_create_anti_targets", {}).get("extra", ""),
    log:
        "references/cnvkit_create_anti_targets/cnvkit_manifest.antitarget.bed.log",
    benchmark:
        repeat(
            "references/cnvkit_create_anti_targets/cnvkit_manifest.antitarget.bed.benchmark.tsv",
            config.get("cnvkit_create_anti_targets", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("cnvkit_create_anti_targets", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("cnvkit_create_anti_targets", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("cnvkit_create_anti_targets", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("cnvkit_create_anti_targets", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("cnvkit_create_anti_targets", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("cnvkit_create_anti_targets", {}).get("time", config["default_resources"]["time"]),
    conda:
        "../envs/cnvkit_panel_of_normal.yaml"
    container:
        config.get("cnvkit_create_anti_targets", {}).get("container", config["default_container"])
    message:
        "{rule}: create anti target bed file for cnvkit from design bed"
    shell:
        "(cnvkit.py antitarget {input.bed} -o {output.bed} {params.extra}) &> {log}"


rule cnvkit_build_normal_reference:
    input:
        bams=lambda wildcards: get_bams(units),
        target="references/cnvkit_create_targets/cnvkit_manifest.target.bed",
        antitarget="references/cnvkit_create_anti_targets/cnvkit_manifest.antitarget.bed",
        ref=config.get("reference", {}).get("fasta", ""),
        mappability=config.get("reference", {}).get("mappability", ""),
    output:
        PoN=temp("references/cnvkit_build_normal_reference/cnvkit.PoN.cnn"),
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
        "references/cnvkit_build_normal_reference/cnvkit.PoN.cnn.log",
    benchmark:
        repeat(
            "references/cnvkit_build_normal_reference/cnvkit.PoN.cnn.benchmark.tsv",
            config.get("cnvkit_build_normal_reference", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("cnvkit_build_normal_reference", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("cnvkit_build_normal_reference", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("cnvkit_build_normal_reference", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("cnvkit_build_normal_reference", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("cnvkit_build_normal_reference", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("cnvkit_build_normal_reference", {}).get("time", config["default_resources"]["time"]),
    conda:
        "../envs/cnvkit_panel_of_normal.yaml"
    container:
        config.get("cnvkit_build_normal_reference", {}).get("container", config["default_container"])
    message:
        "{rule}: build cnvkit normal db"
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
