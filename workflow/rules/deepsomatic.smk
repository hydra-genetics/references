__author__ = "Andrei Guliaev"
__copyright__ = "Copyright 2025, Andrei Guliaev"
__email__ = "andrei.guliaev@scilifelab.uu.se"
__license__ = "GPL-3"


rule deepsomatic_pon:
    input:
        bam=lambda wildcards: get_input_aligned_bam(wildcards, config)[0],
        bai=lambda wildcards: get_input_aligned_bam(wildcards, config)[1],
        ref=config.get("reference", {}).get("fasta", ""),
        bed=config.get("reference", {}).get("design_bed", ""),
    output:
        tmpdir=temp(directory("references/deepsomatic_pon/{sample}.tmp")),
        vcf=temp("references/deepsomatic_pon/{sample}.pon.vcf.gz"),
    params:
        extra=config.get("deepsomatic_pon", {}).get("extra", ""),
        model=config.get("deepsomatic_pon", {}).get("model", ""),
        name=lambda wildcards: f"{wildcards.sample}",
    log:
        dir="references/deepsomatic_pon/{sample}.deepsomatic_pon.dir.log",
        stdout="references/deepsomatic_pon/{sample}.deepsomatic_pon.stdout.log",
    benchmark:
        repeat(
            "references/deepsomatic_pon/{sample}.output.benchmark.tsv",
            config.get("deepsomatic_pon", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("deepsomatic_pon", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("deepsomatic_pon", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("deepsomatic_pon", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("deepsomatic_pon", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("deepsomatic_pon", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("deepsomatic_pon", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("deepsomatic_pon", {}).get("container", config["default_container"])
    message:
        "{rule}: Calling small variants from sequencing data in normal only sample with DeepSomatic from {input.bam}"
    shell:
        """
        run_deepsomatic \
        --model_type={params.model} \
        --ref={input.ref} \
        --reads_tumor={input.bam} \
        --output_vcf={output.vcf} \
        --sample_name_tumor={params.name} \
        --num_shards={resources.threads} \
        --logging_dir={log.dir} \
        --vcf_stats_report=true \
        --intermediate_results_dir {output.tmpdir} \
        --process_somatic=true \
        --regions={input.bed} \
        {params.extra} $> {log.stdout}
        """
