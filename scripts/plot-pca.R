log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library("DESeq2")
library("ggplot2")

# load deseq2 data
dds <- readRDS(snakemake@input[[1]])

# obtain normalized counts
counts <- rlog(dds, blind=FALSE)
svg(snakemake@output[[1]])
plotPCA(counts, intgroup=snakemake@params[["pca_labels"]])
dev.off()


# pcadt <- plotPCA(counts, intgroup=snakemake@params[["pca_labels"]], returnData=T)
# svg(snakemake@output[[2]])
# ggplot(pcadt, aes(PC1, PC2, color = condition)) +
#   geom_point(shape = 16, size = 5, show.legend = TRUE, alpha = 1) +
#   geom_text(aes(label = name), position = nudge)
# dev.off()


# http://www.sthda.com/english/wiki/print.php?id=204
# look at which genes have highest contribution to PCA
# library("FactoMineR")
# library("factoextra")
# res.pca <- PCA(t(assay(counts)), graph = FALSE)
# fviz_contrib(res.pca, choice = "var", axes = 1, top = 50)


# look at pairwise correlation of samples - using DESeq2 norm read counts
# library("PerformanceAnalytics")
# pdf('counts_deseq2tf_pairwise_correlation.pdf')
# chart.Correlation(assay(counts), histogram=TRUE, pch=19)
# dev.off()