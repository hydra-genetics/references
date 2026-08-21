__author__ = "Jonas Almlöf"
__copyright__ = "Copyright 2026, Jonas Almlöf"
__email__ = "jonas.almlof@scilifelab.uu.se"
__license__ = "GPL-3"


rule ichorcna_offtarget_read_counter:
    input:
        bam=lambda wildcards: get_units(units, wildcards)[0].bam,
        bai=lambda wildcards: "%s.bai" % get_units(units, wildcards)[0].bam,
    output:
        wig=temp("references/ichorcna_offtarget_read_counter/{sample}_{type}.wig"),
    params:
        chrs=config.get("ichorcna_offtarget_read_counter", {}).get(
            "chrs",
            "chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,"
            "chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX",
        ),
        window=config.get("ichorcna_offtarget_read_counter", {}).get("window", 100000),
        quality=config.get("ichorcna_offtarget_read_counter", {}).get("quality", 20),
    log:
        "references/ichorcna_offtarget_read_counter/{sample}_{type}.wig.log",
    benchmark:
        repeat(
            "references/ichorcna_offtarget_read_counter/{sample}_{type}.wig.benchmark.tsv",
            config.get("ichorcna_offtarget_read_counter", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("ichorcna_offtarget_read_counter", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("ichorcna_offtarget_read_counter", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("ichorcna_offtarget_read_counter", {}).get(
            "mem_per_cpu", config["default_resources"]["mem_per_cpu"]
        ),
        partition=config.get("ichorcna_offtarget_read_counter", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("ichorcna_offtarget_read_counter", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("ichorcna_offtarget_read_counter", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("ichorcna_offtarget_read_counter", {}).get("container", config["default_container"])
    message:
        "{rule}: Count reads in {params.window}bp bins for {input.bam} with HMMcopy readCounter"
    shell:
        "readCounter {input.bam} -c {params.chrs} -w {params.window} -q {params.quality} > {output.wig} 2> {log}"


rule ichorcna_offtarget_wig_list:
    input:
        wig_files=get_ichorcna_wigs(samples, units),
    output:
        wig_list_file="references/ichorcna_offtarget_wig_list/wig_files.list",
    log:
        "references/ichorcna_offtarget_wig_list/wig_files.list.log",
    benchmark:
        repeat(
            "references/ichorcna_offtarget_wig_list/wig_files.list.benchmark.tsv",
            config.get("ichorcna_offtarget_wig_list", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("ichorcna_offtarget_wig_list", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("ichorcna_offtarget_wig_list", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("ichorcna_offtarget_wig_list", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("ichorcna_offtarget_wig_list", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("ichorcna_offtarget_wig_list", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("ichorcna_offtarget_wig_list", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("ichorcna_offtarget_wig_list", {}).get("container", config["default_container"])
    message:
        "{rule}: make a wig file list for ichorCNA off-target panel of normals"
    shell:
        "(for val in {input.wig_files}; do echo $val >> {output.wig_list_file}; done) &> {log}"


_ichorcna_offtarget_pon_method = config.get("ichorcna_offtarget_panel_of_normals", {}).get("method", "median")


rule ichorcna_offtarget_panel_of_normals:
    input:
        wig_list_file="references/ichorcna_offtarget_wig_list/wig_files.list",
        wig_files=get_ichorcna_wigs(samples, units),
    output:
        pon="references/ichorcna_offtarget_panel_of_normals/ichorcna_offtarget_PoN_%s.rds" % _ichorcna_offtarget_pon_method,
        txt="references/ichorcna_offtarget_panel_of_normals/ichorcna_offtarget_PoN_%s.txt" % _ichorcna_offtarget_pon_method,
    params:
        outfile=lambda wildcards, output: output.pon[: -len("_%s.rds" % _ichorcna_offtarget_pon_method)],
        # gc_wig/map_wig/centromere are params, not input: they're often paths
        # bundled inside the container (e.g. /opt/ichorCNA/inst/extdata/...),
        # which don't exist on the host filesystem Snakemake itself checks
        # against - declaring them as input made dry-run/DAG-building fail
        # with a false "missing input file" for any such path.
        gc_wig=config.get("ichorcna_offtarget_panel_of_normals", {}).get("gc_wig", ""),
        map_wig=config.get("ichorcna_offtarget_panel_of_normals", {}).get("map_wig", ""),
        centromere=config.get("ichorcna_offtarget_panel_of_normals", {}).get("centromere", ""),
        chrs=config.get("ichorcna_offtarget_panel_of_normals", {}).get("chrs", 'c(1:22,"X")'),
        chr_normalize=config.get("ichorcna_offtarget_panel_of_normals", {}).get("chr_normalize", "c(1:22)"),
        # createPanelOfNormals.R never forwards --genomeStyle to its own call to
        # loadReadCountsFromWig(), which defaults to NCBI internally regardless -
        # sample counts get force-converted to NCBI style no matter what is passed
        # here. Defaulting to NCBI keeps chrs/centromere consistent with that.
        genome_style=config.get("ichorcna_offtarget_panel_of_normals", {}).get("genome_style", "NCBI"),
        method=_ichorcna_offtarget_pon_method,
        extra=config.get("ichorcna_offtarget_panel_of_normals", {}).get("extra", ""),
    log:
        "references/ichorcna_offtarget_panel_of_normals/ichorcna_offtarget_PoN.output.log",
    benchmark:
        repeat(
            "references/ichorcna_offtarget_panel_of_normals/ichorcna_offtarget_PoN.output.benchmark.tsv",
            config.get("ichorcna_offtarget_panel_of_normals", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("ichorcna_offtarget_panel_of_normals", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("ichorcna_offtarget_panel_of_normals", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("ichorcna_offtarget_panel_of_normals", {}).get(
            "mem_per_cpu", config["default_resources"]["mem_per_cpu"]
        ),
        partition=config.get("ichorcna_offtarget_panel_of_normals", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("ichorcna_offtarget_panel_of_normals", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("ichorcna_offtarget_panel_of_normals", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("ichorcna_offtarget_panel_of_normals", {}).get("container", config["default_container"])
    message:
        "{rule}: create ichorCNA off-target panel of normals from {input.wig_list_file}"
    shell:
        "(Rscript /opt/ichorCNA/scripts/createPanelOfNormals.R "
        "--filelist {input.wig_list_file} "
        "--gcWig {params.gc_wig} "
        "--mapWig {params.map_wig} "
        "--centromere {params.centromere} "
        "--chrs '{params.chrs}' "
        "--chrNormalize '{params.chr_normalize}' "
        "--genomeStyle {params.genome_style} "
        "--method {params.method} "
        "{params.extra} "
        "--outfile {params.outfile}) &> {log}"
