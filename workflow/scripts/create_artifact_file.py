
import gzip
import statistics
import sys

vcf_files = snakemake.input.vcfs
artifact_panel = open(snakemake.output.artifact_panel, "w")
callers = snakemake.params.callers

FFPE_call_dict = {}
for file_name in vcf_files:
    FFPE_rm_dup_dict = {}
    with gzip.open(file_name, 'rt') as infile:
        file_content = infile.read().split("\n")
        header = True
        for line in file_content:
            if header:
                if line[:6] == "#CHROM":
                    header = False
                continue
            columns = line.strip().split("\t")
            if len(columns) <= 1:
                continue
            chrom = columns[0]
            pos = int(columns[1])
            ref = columns[3]
            alt = columns[4]
            variant_type = "SNV"
            if len(ref) > 1 or len(alt) > 1:
                variant_type = "INDEL"
            INFO = columns[7]
            if INFO[0:3] == "AA=":
                continue
            callers = INFO.split("CALLERS=")[1].split(";")[0].split(",")
            AF = INFO.split(";AF=")
            if len(AF) == 1:
                AF = INFO.split("AF=")
            AF = float(AF.split(";")[0])
            key = chrom + "_" + str(pos) + "_" + variant_type
            if key not in FFPE_call_dict:
                FFPE_call_dict[key] = {}
                for caller in callers:
                    FFPE_call_dict[key][caller] = [0, []]
            if key not in FFPE_rm_dup_dict:
                FFPE_rm_dup_dict[key] = {}
                for caller in callers:
                    FFPE_rm_dup_dict[key][caller] = 0
            for caller in callers:
                if FFPE_rm_dup_dict[key][caller] == 0:
                    FFPE_rm_dup_dict[key][caller] += 1
                    FFPE_call_dict[key][caller][0] += 1
                    FFPE_call_dict[key][caller][1].append(AF)


artifact_panel.write("Chromosome\tpos\tvariant_type")
for caller in callers:
    artifact_panel.write("\t" + caller + "\tmedian_MAF\tsd_MAF")
artifact_panel.write("\n")
for key in FFPE_call_dict:
    chrom = key.split("_")[0]
    pos = key.split("_")[1]
    variant_type = key.split("_")[2]
    artifact_panel.write(chrom + "\t" + pos + "\t" + variant_type)
    for caller in callers:
        FFPE_call_dict[key][caller][1].sort()
        median_af = statistics.median(FFPE_call_dict[key][caller][1])
        sd_af = 1000
        nr_obs = len(FFPE_call_dict[key][caller][1])
        if nr_obs >= 4:
            '''This is the sample variance s² with Bessel’s correction, also known as variance with N-1 degrees of freedom.
            Provided that the data points are representative (e.g. independent and identically distributed),
            the result should be an unbiased estimate of the true population variance.'''
            sd_af = statistics.stdev(FFPE_call_dict[key][caller][1])
        artifact_panel.write("\t" + str(FFPE_call_dict[key][caller]) + "\t" + str(median_af) + "\t" + str(sd_af))
    artifact_panel.write("\n")
