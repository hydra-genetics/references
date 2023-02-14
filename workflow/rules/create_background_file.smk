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
        "references/create_background_file/background_panel.tsv.log",
    benchmark:
        repeat(
            "references/create_background_file/background_panel.tsv.benchmark.tsv",
            config.get("create_background_file", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("create_background_file", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("create_background_file", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("create_background_file", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("create_background_file", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("create_background_file", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("create_background_file", {}).get("time", config["default_resources"]["time"]),
    conda:
        "../envs/create_background_file.yaml"
    container:
        config.get("create_background_file", {}).get("container", config["default_container"])
    message:
        "{rule}: create background PoN"
    script:
        "../scripts/create_background_file.py"
