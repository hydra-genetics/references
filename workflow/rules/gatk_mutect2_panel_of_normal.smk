__author__ = "Jessika Nordin"
__copyright__ = "Copyright 2025, Jessika Nordin"
__email__ = "jessika.nordin@scilifelab.uu.se"
__license__ = "GPL-3"


rule create_somatic_pon:
    input:
        dir="reference/genomicsdbimport/pon_db",  #### Is this needed?
        gendb="gendb://reference/genomicsdbimport/pon_db",
    output:
        vcf="reference/create_somatic_pon/pon.vcf.gz",
    params:
        extra=config.get("create_somatic_pon", {}).get("extra", ""),
    log:
        "reference/create_somatic_pon/pon.vcf.gz.log",
    benchmark:
        repeat(
            "reference/create_somatic_pon/pon.vcf.gz.benchmark.tsv",
            config.get("create_somatic_pon", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("create_somatic_pon", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("create_somatic_pon", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("create_somatic_pon", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("create_somatic_pon", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("create_somatic_pon", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("create_somatic_pon", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("create_somatic_pon", {}).get("container", config["default_container"])
    message:
        "{rule}: generate a normal panel from mutect2 using gvcfs"
    shell:
        "(gatk CreateSomaticPanelOfNormals "
        "-V {input.gendb} "
        "-O {output.vcf} "
        "{params.extra}) &> {log}"


rule genomicsdbimport:
    input:
        fasta=config.get("reference", {}).get("fasta", ""),
        gvcfs=get_gvcfs(units),
        interval=config.get("reference", {}).get("design_bed", ""),
        normal_list=config.get("reference", {}).get("normals_list", ""),
    output:
        gendb="gendb://reference/genomicsdbimport/pon_db",
        path="reference/genomicsdbimport/pon_db",
    params:
        extra=config.get("genomicsdbimport", {}).get("extra", ""),
    log:
        "reference/genomicsdbimport/genomicsdbimport.vcf.gz.log",
    benchmark:
        repeat(
            "reference/genomicsdbimport/genomicsdbimport.vcf.gz.benchmark.tsv",
            config.get("genomicsdbimport", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("genomicsdbimport", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("genomicsdbimport", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("genomicsdbimport", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("genomicsdbimport", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("genomicsdbimport", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("genomicsdbimport", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("genomicsdbimport", {}).get("container", config["default_container"])
    message:
        "{rule}: generate a normal panel from mutect2 using gvcfs"
    shell:
        "(gatk GenomicsDBImport "
        "--genomicsdb-workspace-path {output.path} "
        "--L {input.interval} "
        "-R {input.fasta} "
        "--reader-threads {resources.threads} "
        "--sample-name-map {input.normal_list} "
        "{params.extra}) &> {log}"
