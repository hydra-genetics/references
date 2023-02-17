__author__ = "Jonas Almlöf, Patrik Smeds"
__copyright__ = "Copyright 2021, Patrik Smeds, Jonas Almlöf"
__email__ = "jonas.almlöf@scilifelab.uu.se, patrik.smeds@scilifelab.uu.se"
__license__ = "GPL3"


rule create_artifact_file:
    input:
        vcfs=get_vcfs(units),
    output:
        artifact_panel=temp("references/create_artifact_file/artifact_panel.tsv"),
    params:
        callers=config.get("create_artifact_file", {}).get("callers", ["vardict", "gatk_mutect2"]),
    log:
        "references/create_artifact_file/artifact_panel.tsv.log",
    benchmark:
        repeat(
            "references/create_artifact_file/artifact_panel.tsv.benchmark.tsv",
            config.get("create_artifact_file", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("create_artifact_file", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("create_artifact_file", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("create_artifact_file", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("create_artifact_file", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("create_artifact_file", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("create_artifact_file", {}).get("time", config["default_resources"]["time"]),
    conda:
        "../envs/create_artifact_file.yaml"
    container:
        config.get("cnvkit_create_targets", {}).get("container", config["default_container"])
    message:
        "{rule}: create artifact PoN"
    script:
        "../scripts/create_artifact_file.py"
