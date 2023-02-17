library(ExomeDepth)

bam_df <- read.delim(snakemake@input[["bam_list_file"]], header = FALSE)
bam_list <- bam_df$V1
exons_bed <- read.csv(snakemake@params[["bed"]], header = FALSE, sep = "\t")

# Create counts dataframe for all BAMs
my_counts <- getBamCounts(bed.frame = exons_bed,
                          bam.files = bam_list,
                          include.chr = FALSE)

# prepare the main matrix of read count data
refcount_mat <- as.matrix(
    my_counts[, grep(names(my_counts), pattern = ".bam$")])
save(refcount_mat, file = snakemake@output[["reference"]])
