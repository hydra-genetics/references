import gzip
import statistics

vcf_files = snakemake.input.vcfs
artifact_panel_path = snakemake.output.artifact_panel
callers_used = snakemake.params.callers

FFPE_call_dict = {}


def get_info_field(info_str, key):
    """Safely extracts a value from the INFO field string."""
    for field in info_str.split(";"):
        if field.startswith(f"{key}="):
            return field.split("=", 1)[1]
    return None


for file_name in vcf_files:
    FFPE_rm_dup_dict = {}
    with gzip.open(file_name, "rt") as infile:
        for line in infile:
            if line.startswith("#"):
                continue

            columns = line.strip().split("\t")
            if len(columns) < 8:
                continue

            chrom = columns[0]
            pos = int(columns[1])
            ref = columns[3]
            alt = columns[4]
            variant_type = "SNV" if len(ref) == 1 and len(alt) == 1 else "INDEL"

            info = columns[7]
            if get_info_field(info, "AA") is not None:
                continue

            callers_str = get_info_field(info, "CALLERS")
            if not callers_str:
                continue
            callers = callers_str.split(",")

            af_str = get_info_field(info, "AF")
            if not af_str:
                continue
            try:
                af = float(af_str.split(",")[0])  # Handle potential multi-allelic AF
            except (ValueError, IndexError):
                continue

            key = (chrom, pos, variant_type)

            if key not in FFPE_call_dict:
                FFPE_call_dict[key] = {c: [0, []] for c in callers}
            else:
                for caller in callers:
                    if caller not in FFPE_call_dict[key]:
                        FFPE_call_dict[key][caller] = [0, []]

            if key not in FFPE_rm_dup_dict:
                FFPE_rm_dup_dict[key] = {c: 0 for c in callers}
            else:
                for caller in callers:
                    if caller not in FFPE_rm_dup_dict[key]:
                        FFPE_rm_dup_dict[key][caller] = 0

            for caller in callers:
                if FFPE_rm_dup_dict[key][caller] == 0:
                    FFPE_rm_dup_dict[key][caller] = 1
                    FFPE_call_dict[key][caller][0] += 1
                    FFPE_call_dict[key][caller][1].append(af)
                else:
                    # If multiallelic, use the variant with the highest AF
                    if af > FFPE_call_dict[key][caller][1][-1]:
                        FFPE_call_dict[key][caller][1][-1] = af

with open(artifact_panel_path, "w") as artifact_panel:
    header = ["Chromosome", "pos", "variant_type"]
    for caller in callers_used:
        header.extend([caller, "median_MAF", "sd_MAF"])
    artifact_panel.write("\t".join(header) + "\n")

    for key, callers_data in FFPE_call_dict.items():
        chrom, pos, variant_type = key
        row = [str(chrom), str(pos), str(variant_type)]
        for caller in callers_used:
            if caller not in callers_data:
                row.extend(["0", "0", "1000"])
            else:
                obs_count, af_list = callers_data[caller]
                af_list.sort()
                median_af = statistics.median(af_list)
                sd_af = 1000
                if len(af_list) >= 4:
                    sd_af = statistics.stdev(af_list)
                row.extend([str(obs_count), str(median_af), str(sd_af)])
        artifact_panel.write("\t".join(row) + "\n")
