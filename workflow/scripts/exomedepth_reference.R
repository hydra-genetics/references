library(ExomeDepth)

bam_df <- read.delim(snakemake@input[["bam_list_file"]], header = FALSE)
bam_list <- bam_df$V1
exons_bed <- read.delim(snakemake@params[["bed"]], header = FALSE, sep = "\t")

# Create counts dataframe for all BAMs
refcount_df <- getBamCounts(bed.frame = exons_bed,
                          bam.files = bam_list,
                          include.chr = FALSE)

save(refcount_df, file = snakemake@output[["reference"]])
