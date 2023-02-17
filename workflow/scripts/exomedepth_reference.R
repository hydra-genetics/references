library(ExomeDepth)

bam_list <- read.delim(snakemake@input[["bam_list_file"]], header = FALSE)
exons_bed <- read.csv(snakemake@input[["bed"]], header = FALSE, sep = "\t")

# Create counts dataframe for all BAMs
my_counts <- getBamCounts(bed.frame = exons_bed,
                          bam.files = bam_list,
                          include.chr = FALSE)

# prepare the main matrix of read count data
refcount_mat <- as.matrix(
    my.counts[, grep(names(my.counts), pattern = ".bam$")])
save(refcount_mat, file = snakemake@output[["reference"]])
