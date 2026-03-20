import sys

bams = snakemake.input.bams
PoN_list = open(snakemake.output.PoN_list, "w")

for bam in bams:
    sample = bam.split("/")[-1]
    PoN_list.write(sample + "\t" + bam + "\n")
PoN_list.close()
