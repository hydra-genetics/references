__author__ = "Jonas Almlöf"
__copyright__ = "Copyright 2022, Jonas Almlöf"
__email__ = "jonas.almlöf@scilifelab.uu.se"
__license__ = "GPL3"


rule purecn_interval_file:
    input:
        ref_fasta=config.get("reference", {}).get("fasta", ""),
        design_bed=config.get("reference", {}).get("design_bed", ""),
    output:
        intervals_file=temp("references/purecn_interval_file/targets_intervals.txt"),
        optimized_bed=temp("references/purecn_interval_file/targets_optimized.bed"),
    params:
        genome=config.get("purecn_interval_file", {}).get("genome", "hg19"),
        average_off_target_width=config.get("purecn_interval_file", {}).get("average_off_target_width", "25000"),
    log:
        "references/purecn_interval_file/targets.log",
    benchmark:
        repeat(
            "references/purecn_interval_file/targets.benchmark.tsv",
            config.get("purecn_interval_file", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("purecn_interval_file", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("purecn_interval_file", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("purecn_interval_file", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("purecn_interval_file", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("purecn_interval_file", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("purecn_interval_file", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("purecn_interval_file", {}).get("container", config["default_container"])
    conda:
        "../envs/purecn.yaml"
    message:
        "{rule}: make an interval file for purecn"
    shell:
        "(Rscript $PURECN/IntervalFile.R "
        "--fasta {input.ref_fasta} "
        "--in-file {input.design_bed} "
        "--out-file {output.intervals_file} "
        "--export {output.optimized_bed} "
        "--genome {params.genome} "
        "--average-off-target-width {params.average_off_target_width} "
        "--off-target) &> {log}"


rule purecn_bam_list:
    input:
        bam_list=get_bams(units),
    output:
        bam_list_file="references/purecn_bam_list/bam_files.list",
    log:
        "references/purecn_bam_list/bam_files.list.log",
    benchmark:
        repeat(
            "references/purecn_bam_list/bam_files.list.benchmark.tsv",
            config.get("purecn_bam_list", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("purecn_bam_list", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("purecn_bam_list", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("purecn_bam_list", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("purecn_bam_list", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("purecn_bam_list", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("purecn_bam_list", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("purecn_bam_list", {}).get("container", config["default_container"])
    conda:
        "../envs/purecn.yaml"
    message:
        "{rule}: make a bam list file for purecn coverage"
    shell:
        "(for val in {input.bam_list}; do echo $val >> {output.bam_list_file}; done) &> {log}"


rule purecn_coverage:
    input:
        bam_list_file="references/purecn_bam_list/bam_files.list",
    output:
        coverage_list=get_coverage_files(samples, units),
    params:
        intervals=config.get("purecn_coverage", {}).get("intervals", ""),
        extra=config.get("purecn_coverage", {}).get("extra", ""),
    log:
        "references/purecn_coverage/purecn_coverage.output.log",
    benchmark:
        repeat(
            "references/purecn_coverage/purecn_coverage.output.benchmark.tsv",
            config.get("purecn_coverage", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("purecn_coverage", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("purecn_coverage", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("purecn_coverage", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("purecn_coverage", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("purecn_coverage", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("purecn_coverage", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("purecn_coverage", {}).get("container", config["default_container"])
    conda:
        "../envs/purecn.yaml"
    message:
        "{rule}: calculate coverage for all samples in {input}"
    shell:
        "(Rscript $PURECN/Coverage.R "
        "--out-dir=references/purecn_coverage "
        "--bam={input.bam_list_file} "
        "--intervals={params.intervals} "
        "--cores={threads} "
        "{params.extra}) &> {log}"


rule purecn_coverage_list:
    input:
        coverage_list=get_coverage_files(samples, units),
    output:
        coverage_list_file="references/purecn_coverage_list/coverage_files.list",
    log:
        "references/purecn_coverage_list/coverage_files.list.log",
    benchmark:
        repeat(
            "references/purecn_coverage_list/coverage_files.list.benchmark.tsv",
            config.get("purecn_coverage_list", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("purecn_coverage_list", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("purecn_coverage_list", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("purecn_coverage_list", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("purecn_coverage_list", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("purecn_coverage_list", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("purecn_coverage_list", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("purecn_coverage_list", {}).get("container", config["default_container"])
    conda:
        "../envs/purecn.yaml"
    message:
        "{rule}: make a coverage list file for purecn normalDB"
    shell:
        "(for val in {input.coverage_list}; do echo $val >> {output.coverage_list_file}; done) &> {log}"


rule bcftools_merge:
    input:
        vcfs=get_vcfs(units),
        vcfs_tabix=expand("{dataset}.{ext}", dataset=get_vcfs(units), ext=["tbi"]),
    output:
        normal_vcf="references/bcftools_merge/normal_db.vcf.gz",
    params:
        output_type=config.get("bcftools_merge", {}).get("output_type", "z"),
        info_rules=config.get("bcftools_merge", {}).get("output_type", "-"),
        extra=config.get("bcftools_merge", {}).get("extra", ""),
    log:
        "references/bcftools_merge/normal_db.vcf.log",
    benchmark:
        repeat(
            "references/bcftools_merge/normal_db.vcf.benchmark.tsv",
            config.get("bcftools_merge", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("bcftools_merge", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("bcftools_merge", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("bcftools_merge", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("bcftools_merge", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("bcftools_merge", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("bcftools_merge", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("bcftools_merge", {}).get("container", config["default_container"])
    conda:
        "../envs/purecn.yaml"
    message:
        "{rule}: merge normal sample vcf files into one big vcf file"
    shell:
        "(bcftools merge "
        "{input.vcfs} "
        "--output-type {params.output_type} "
        "--info-rules {params.info_rules} "
        "{params.extra} > {output.normal_vcf}) &> {log}"


rule purecn_normal_db:
    input:
        coverage_list_file="references/purecn_coverage_list/coverage_files.list",
        normal_vcf="references/bcftools_merge/normal_db.vcf.gz",
        normal_vcf_tbi="references/bcftools_merge/normal_db.vcf.gz.tbi",
    output:
        normal_db="references/purecn_normal_db/output/normalDB.rds",
        mapping_bias="references/purecn_normal_db/output/mapping_bias.rds",
        out_dir=directory("references/purecn_normal_db/output/"),
    params:
        genome=config.get("purecn_normal_db", {}).get("genome", "hg19"),
    log:
        "references/purecn_normal_db/normal_db.output.log",
    benchmark:
        repeat(
            "references/purecn_normal_db/normal_db.output.benchmark.tsv",
            config.get("purecn_normal_db", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("purecn_normal_db", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("purecn_normal_db", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("purecn_normal_db", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("purecn_normal_db", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("purecn_normal_db", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("purecn_normal_db", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("purecn_normal_db", {}).get("container", config["default_container"])
    conda:
        "../envs/purecn.yaml"
    message:
        "{rule}: calculate normal DB for purecn"
    shell:
        "(Rscript $PURECN/NormalDB.R "
        "--out-dir {output.out_dir} "
        "--coverage-files {input.coverage_list_file} "
        "--normal-panel {input.normal_vcf} "
        "--genome {params.genome}) &> {log}"
