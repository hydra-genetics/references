__author__ = "Jonas Almlöf, Patrik Smeds"
__copyright__ = "Copyright 2021, Patrik Smeds, Jonas Almlöf"
__email__ = "jonas.almlöf@scilifelab.uu.se, patrik.smeds@scilifelab.uu.se"
__license__ = "GPL3"


rule bed_to_interval_list:
    input:
        bed=config.get("reference", {}).get("design_bed", ""),
        ref=config.get("reference", {}).get("fasta", ""),
        ref_dict=config.get("reference", {}).get("dict", ""),
    output:
        temp("references/bed_to_interval_list/%s.interval_list" % config["reference"]["design_bed"].split("/")[-1]),
    log:
        "references/bed_to_interval_list/bed_to_interval_list.output.log",
    benchmark:
        repeat(
            "references/bed_to_interval_list/bed_to_interval_list.output.benchmark.tsv",
            config.get("bed_to_interval_list", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("bed_to_interval_list", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("bed_to_interval_list", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("bed_to_interval_list", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("bed_to_interval_list", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("bed_to_interval_list", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("bed_to_interval_list", {}).get("time", config["default_resources"]["time"]),
    conda:
        "../envs/gatk_cnv_panel_of_normal.yaml"
    container:
        config.get("bed_to_interval_list", {}).get("container", config["default_container"])
    message:
        "{rule}: convert design bed file to interval file format"
    shell:
        "(gatk BedToIntervalList  -I {input.bed} -O {output} -SD {input.ref} ) &> {log} "


rule preprocess_intervals:
    input:
        ref=config.get("reference", {}).get("fasta", ""),
        intervalList="references/bed_to_interval_list/%s.interval_list" 
        % config.get("reference", {}).get("design_bed", "").split("/")[-1],
    output:
        temp(
            "references/preprocess_intervals/%s.preprocessed.interval_list"
            % config.get("reference", {}).get("design_bed", "").split("/")[-1]
        ),
    params:
        bin_length=config.get("preprocess_intervals", {}).get("bin_length", "0"),  # WGS 1000 Exomes/target: 0
        padding=config.get("preprocess_intervals", {}).get("padding", "250"),  # WGS 0 Exomes/target: 250
        extra=config.get("preprocess_intervals", {}).get("extra", "-imr OVERLAPPING_ONLY"),
    log:
        "references/preprocess_intervals/design.preprocessed.interval_list.log",
    benchmark:
        repeat(
            "references/preprocess_intervals/background_panel.tsv.benchmark.tsv",
            config.get("preprocess_intervals", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("preprocess_intervals", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("preprocess_intervals", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("preprocess_intervals", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("preprocess_intervals", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("preprocess_intervals", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("preprocess_intervals", {}).get("time", config["default_resources"]["time"]),
    conda:
        "../envs/gatk_cnv_panel_of_normal.yaml"
    container:
        config.get("preprocess_intervals", {}).get("container", config["default_container"])
    message:
        "{rule}: preprocess interval file"
    shell:
        "(gatk --java-options '-Xmx4g' PreprocessIntervals -L {input.intervalList} -R {input.ref} "
        "--bin-length {params.bin_length} --padding {params.padding} {params.extra} "
        "-O {output}) &> {log}"


rule collect_read_counts:
    input:
        bam=lambda wildcards: get_units(units, wildcards)[0].bam,
        bai=lambda wildcards: "%s.bai" % get_units(units, wildcards)[0].bam,
        interval="references/preprocess_intervals/%s.preprocessed.interval_list"
        % config.get("reference", {}).get("design_bed", "").split("/")[-1],
    output:
        temp("references/collect_read_counts/{sample}_{type}.counts.hdf5"),
    params:
        extra=config.get("collect_read_counts", {}).get("extra", "-imr OVERLAPPING_ONLY"),
    log:
        "references/collect_read_counts/{sample}_{type}.counts.hdf5.log",
    benchmark:
        repeat(
            "references/collect_read_counts/{sample}_{type}.counts.hdf5.benchmark.tsv",
            config.get("collect_read_counts", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("collect_read_counts", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("collect_read_counts", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("collect_read_counts", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("collect_read_counts", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("collect_read_counts", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("collect_read_counts", {}).get("time", config["default_resources"]["time"]),
    conda:
        "../envs/gatk_cnv_panel_of_normal.yaml"
    container:
        config.get("collect_read_counts", {}).get("container", config["default_container"])
    message:
        "{rule}: collect read counts for each sample"
    shell:
        "(gatk --java-options '-Xmx4g' CollectReadCounts -I {input.bam} -L {input.interval} "
        "{params.extra} -O {output}) &> {log}"


rule create_read_count_panel_of_normals:
    input:
        bams=[
            "references/collect_read_counts/%s_%s.counts.hdf5" % (sample, t)
            for sample in get_samples(samples)
            for t in get_unit_types(units, sample)
        ],
    output:
        temp("references/create_read_count_panel_of_normals/gatk_cnv_panel_of_normal.hdf5"),
    params:
        extra=config.get("create_read_count_panel_of_normals", {}).get("extra", ""),
        input=lambda wildcards, input: " -I ".join(input.bams),
    log:
        "references/create_read_count_panel_of_normals/gatk_cnv_panel_of_normal.hdf5.log",
    benchmark:
        repeat(
            "references/create_read_count_panel_of_normals/gatk_cnv_panel_of_normal.hdf5.benchmark.tsv",
            config.get("create_read_count_panel_of_normals", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("create_read_count_panel_of_normals", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("create_read_count_panel_of_normals", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("create_read_count_panel_of_normals", {}).get(
            "mem_per_cpu", config["default_resources"]["mem_per_cpu"]
        ),
        partition=config.get("create_read_count_panel_of_normals", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("create_read_count_panel_of_normals", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("create_read_count_panel_of_normals", {}).get("time", config["default_resources"]["time"]),
    conda:
        "../envs/gatk_cnv_panel_of_normal.yaml"
    container:
        config.get("create_read_count_panel_of_normals", {}).get("container", config["default_container"])
    message:
        "{rule}: create GATK CNV PoN"
    shell:
        "(gatk --java-options '-Xmx4g' CreateReadCountPanelOfNormals -I {params.input} "
        "{params.extra} "
        "-O {output}) &> {log}"
