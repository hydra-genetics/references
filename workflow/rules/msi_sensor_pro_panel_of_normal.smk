__author__ = "Jonas Almlöf, Patrik Smeds"
__copyright__ = "Copyright 2021, Patrik Smeds, Jonas Almlöf"
__email__ = "jonas.almlöf@scilifelab.uu.se, patrik.smeds@scilifelab.uu.se"
__license__ = "GPL3"


rule msisensor_pro_scan:
    input:
        ref=config.get("reference", {}).get("fasta", ""),
    output:
        PoN_list=temp("references/msisensor_pro_scan/Msisensor_pro_reference.list"),
    log:
        "references/msisensor_pro_scan/Msisensor_pro_reference.list.log",
    benchmark:
        repeat(
            "references/msisensor_pro_scan/Msisensor_pro_reference.list.benchmark.tsv",
            config.get("msisensor_pro_scan", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("msisensor_pro_scan", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("msisensor_pro_scan", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("msisensor_pro_scan", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("msisensor_pro_scan", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("msisensor_pro_scan", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("msisensor_pro_scan", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("msisensor_pro_scan", {}).get("container", config["default_container"])
    message:
        "{rule}: scan reference genome"
    shell:
        "(msisensor-pro scan -d {input.ref} -o {output.PoN_list}) &> {log}"


rule msisensor_pro_input_file:
    input:
        bams=get_bams(units),
    output:
        PoN_list=temp("references/msisensor_pro_input_file/configure.txt"),
    log:
        "references/msisensor_pro_input_file/configure.txt.log",
    benchmark:
        repeat(
            "references/msisensor_pro_input_file/configure.txt.benchmark.tsv",
            config.get("msisensor_pro_input_file", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("msisensor_pro_input_file", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("msisensor_pro_input_file", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("msisensor_pro_input_file", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("msisensor_pro_input_file", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("msisensor_pro_input_file", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("msisensor_pro_input_file", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("msisensor_pro_input_file", {}).get("container", config["default_container"])
    message:
        "{rule}: make msisensor run configuration file"
    script:
        "../scripts/msisensor_pro_input_file.py"


rule msisensor_pro_baseline:
    input:
        bam_conf="references/msisensor_pro_input_file/configure.txt",
        PoN_list="references/msisensor_pro_scan/Msisensor_pro_reference.list",
    output:
        PoN_list=temp("references/msisensor_pro_baseline/Msisensor_pro_reference.list_baseline"),
    params:
        extra=config.get("collect_read_counts", {}).get("extra", "-c 50"),  # -c = minimal coverage, WXS: 20; WGS: 15
        outdir=lambda wildcards, output: os.path.dirname(os.path.abspath(output.PoN_list)),
    log:
        "references/msisensor_pro_baseline/Msisensor_pro_reference.list_baseline.log",
    benchmark:
        repeat(
            "references/msisensor_pro_baseline/Msisensor_pro_reference.list_baseline.benchmark.tsv",
            config.get("msisensor_pro_baseline", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("msisensor_pro_baseline", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("msisensor_pro_baseline", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("msisensor_pro_baseline", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("msisensor_pro_baseline", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("msisensor_pro_baseline", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("msisensor_pro_baseline", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("msisensor_pro_baseline", {}).get("container", config["default_container"])
    message:
        "{rule}: create MSI-sensor pro PoN baseline"
    shell:
        "(msisensor-pro baseline {params.extra} -d {input.PoN_list} -i {input.bam_conf} -o {params.outdir}) &> {log}"
