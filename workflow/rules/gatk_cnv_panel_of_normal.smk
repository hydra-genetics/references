# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

__author__ = "Jonas Almlöf, Patrik Smeds"
__copyright__ = "Copyright 2021, Patrik Smeds, Jonas Almlöf"
__email__ = "jonas.almlöf@scilifelab.uu.se, patrik.smeds@scilifelab.uu.se"
__license__ = "GPL3"


rule bed_to_interval_list:
    input:
        bed=config["reference"]["design_bedfile"],
        ref=config["reference"]["fasta"],
        ref_dict=config["reference"]["dict"],
    output:
        temp("references/bed_to_interval_list/%s.interval_list" % config["reference"]["design_bedfile"].split("/")[-1]),
    log:
        "references/bed_to_interval_list/bedToIntervalList.log",
    conda:
        "../envs/gatk_cnv_panel_of_normal.yaml"
    container:
        config.get("bed_to_interval_list", {}).get("container", config["default_container"])
    shell:
        "(gatk BedToIntervalList  -I {input.bed} -O {output} -SD {input.ref} ) &> {log} "


rule preprocess_intervals:
    input:
        ref=config["reference"]["fasta"],
        intervalList="references/bed_to_interval_list/%s.interval_list" % config["reference"]["design_bedfile"].split("/")[-1],
    output:
        "references/preprocess_intervals/%s.preprocessed.interval_list" % config["reference"]["design_bedfile"].split("/")[-1],
    params:
        bin_length=config.get("preprocess_intervals", {}).get("bin_length", "0"),  # WGS 1000 Exomes/target: 0
        padding=config.get("preprocess_intervals", {}).get("padding", "250"),  # WGS 0 Exomes/target: 250
        extra=config.get("preprocess_intervals", {}).get("extra", "-imr OVERLAPPING_ONLY"),
    log:
        "references/preprocess_intervals/design.preprocessed.interval_list.log",
    conda:
        "../envs/gatk_cnv_panel_of_normal.yaml"
    container:
        config.get("preprocess_intervals", {}).get("container", config["default_container"])
    shell:
        "(gatk --java-options '-Xmx4g' PreprocessIntervals -L {input.intervalList} -R {input.ref} "
        "--bin-length {params.bin_length} --padding {params.padding} {params.extra} "
        "-O {output}) &> {log}"


rule collect_read_counts:
    input:
        bam=lambda wildcards: get_units(units, wildcards)[0].bam,
        bai=lambda wildcards: "%s.bai" % get_units2(units, wildcards)[0].bam,
        interval="references/preprocess_intervals/%s.preprocessed.interval_list"
        % config["reference"]["design_bedfile"].split("/")[-1],
    output:
        temp("references/collect_read_counts/{sample}_{type}.counts.hdf5"),
    params:
        extra=config.get("collect_read_counts", {}).get("extra", "-imr OVERLAPPING_ONLY"),
    log:
        "references/collect_read_counts/{sample}_{type}.counts.hdf5.log",
    conda:
        "../envs/gatk_cnv_panel_of_normal.yaml"
    container:
        config.get("collect_read_counts", {}).get("container", config["default_container"])
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
        "references/create_read_count_panel_of_normals/gatk_cnv_panel_of_normal.hdf5",
    params:
        extra=config.get("create_read_count_panel_of_normals", {}).get("extra", ""),
        input=lambda wildcards, input: " -I ".join(input.bams),
    log:
        "references/create_read_count_panel_of_normals/GATK/gatk_cnv_panel_of_normal.hdf5.log",
    conda:
        "../envs/gatk_cnv_panel_of_normal.yaml"
    container:
        config.get("create_read_count_panel_of_normals", {}).get("container", config["default_container"])
    shell:
        "(gatk --java-options '-Xmx4g' CreateReadCountPanelOfNormals -I {params.input} "
        "{params.extra} "
        "-O {output}) &> {log}"
