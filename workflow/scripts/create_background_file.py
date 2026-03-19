import gzip
import statistics


gvcf_files = snakemake.input.gvcfs
background_file_path = snakemake.output.background_file
min_dp = snakemake.params.min_dp
max_af = snakemake.params.max_af

background_dict = {}

for file_name in gvcf_files:
    with gzip.open(file_name, "rt") as infile:
        for line in infile:
            if line.startswith("#"):
                continue

            columns = line.strip().split("\t")
            if len(columns) < 10:
                continue

            # Handle potential chromosome prefixes
            chrom = columns[0]
            if chrom.startswith("chr"):
                chrom = chrom[3:]

            pos = columns[1]
            key = (chrom, pos)

            formats = columns[8].split(":")
            data = columns[9].split(":")

            try:
                ad_idx = formats.index("AD")
                ad_info = data[ad_idx].split(",")

                ref_ad = int(ad_info[0])
                alt_ad = sum(int(ad) for ad in ad_info[1:])
                dp = ref_ad + alt_ad

                if dp <= min_dp:
                    continue

                alt_af = alt_ad / float(dp)

                # Check if it's potentially a homozygous variant being filtered
                if alt_af > 1 - max_af:
                    alt_af = 1 - alt_af

                if alt_af > max_af:
                    continue

                if key in background_dict:
                    background_dict[key].append(alt_af)
                else:
                    background_dict[key] = [alt_af]
            except (ValueError, IndexError):
                # Skip if AD field or index is missing/malformed
                continue

with open(background_file_path, "w") as background_file:
    background_file.write("Chr\tPos\tMedian\tSD\tNrObs\n")
    for key, af_list in background_dict.items():
        if len(af_list) < 4:
            continue

        median_background = statistics.median(af_list)
        stdev_background = statistics.stdev(af_list)

        chrom, pos = key
        background_file.write(f"{chrom}\t{pos}\t{median_background}\t{stdev_background}\t{len(af_list)}\n")
