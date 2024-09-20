__author__ = "Jonas Almlöf"
__copyright__ = "Copyright 2024, Jonas Almlöf"
__email__ = "jonas.almlof@scilifelab.uu.se"
__license__ = "GPL-3"


rule jumble_count:
    input:
        bam=lambda wildcards: get_units(units, wildcards)[0].bam,
        bai=lambda wildcards: "%s.bai" % get_units(units, wildcards)[0].bam,
    output:
        counts=temp("alignment/jumble_count/{sample}_{type}.bam.counts.RDS"),
    params:
        bed=config.get("reference", {}).get("design_bed", ""),
    log:
        "references/jumble_count/{sample}_{type}.output.log",
    benchmark:
        repeat(
            "references/jumble_count/{sample}_{type}.output.benchmark.tsv",
            config.get("jumble_count", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("jumble_count", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("jumble_count", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("jumble_count", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("jumble_count", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("jumble_count", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("jumble_count", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("jumble_count", {}).get("container", config["default_container"])
    message:
        "{rule}: generate counts in for supplied bedfile in {input.bam}"
    shell:
        "(Rscript /Jumble/jumble-count.R "
        "-t {params.bed} "
        "-b {input.bam}) &> {log}"


rule jumble_reference:
    input:
        count_files=get_counts(samples, units),
    output:
        PoN="references/jumble_reference/%s.reference.RDS" % config.get("reference", {}).get("design_bed", "").split("/")[-1],
    params:
        annotation=config.get("jumble_reference", {}).get("annotation", ""),
        bed=config.get("reference", {}).get("design_bed", ""),
        input_dir=lambda wildcards, input: os.path.dirname(input[0]),
        output_dir="alignment/jumble_reference/",
    log:
        "references/jumble_reference/%s.reference.RDS.output.log"
        % config.get("reference", {}).get("design_bed", "").split("/")[-1],
    benchmark:
        repeat(
            f"references/jumble_reference/%s.reference.RDS.output.benchmark.tsv"
            % config.get("reference", {}).get("design_bed", "").split("/")[-1],
            config.get("jumble_reference", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("jumble_reference", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("jumble_reference", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("jumble_reference", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("jumble_reference", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("jumble_reference", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("jumble_reference", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("jumble_reference", {}).get("container", config["default_container"])
    message:
        "{rule}: from jumble count files make a PoN {output.PoN}"
    shell:
        "(Rscript /Jumble/jumble-reference.R "
        "-i {params.input_dir} "
        "-a {params.annotation} "
        "-o {params.output_dir} "
        "-c {threads}) &> {log}"
