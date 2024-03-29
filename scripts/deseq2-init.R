log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library("DESeq2")

parallel <- FALSE
if (snakemake@threads > 1) {
    library("BiocParallel")
    # setup parallelization
    register(MulticoreParam(snakemake@threads))
    parallel <- TRUE
}

# colData and countData must have the same sample order, but this is ensured
# by the way we create the count matrix
cts <- read.table(snakemake@input[["counts"]], header=TRUE, row.names="gene")
coldata <- read.table(snakemake@input[["samples"]], header=TRUE, row.names="sample")

colnames(cts) <- rownames(coldata)
dds <- DESeqDataSetFromMatrix(countData=cts,
                              colData=coldata["condition"],
                              design=~ condition)

# remove uninformative columns
dds <- dds[ rowSums(counts(dds)) > 1, ]
# TODO optionally allow to collapse technical replicates
dds <- DESeq(dds, parallel=parallel)

saveRDS(dds, file=snakemake@output[[1]])
