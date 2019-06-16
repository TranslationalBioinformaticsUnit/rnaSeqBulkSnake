log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library("edgeR")
library("Biobase")
library("SummarizedExperiment")

# inspired by 
# https://github.com/mikelove/deseq2_or_edger/blob/master/runScripts.R

# recommended filter: Gordon Smyth
# https://support.bioconductor.org/p/85511/#85514
cpmFilter <- function(e, cpm=TRUE) {
  if (cpm) {
    L <- min(colSums(exprs(e))/1e6)
    dgel <- DGEList(exprs(e))
    # here I use 3 even for the heldout set (n=7/8),
    # because otherwise the filtering is too strict
    # to the detriment of edgeR and limma-voom performance
    return(rowSums(cpm(dgel) > 10/L) >= 3)
  } else {
    return(rowSums(exprs(e)) > 0)
  }
}

# Genewise Negative Binomial Generalized Linear Models
runEdgeR <- function(e, coldt) {
  padj <- rep(NA, nrow(e))
  pval <- rep(NA, nrow(e))
  logFc <- rep(NA, nrow(e))

  keep <- cpmFilter(e, cpm=TRUE)
  gene.names <- rownames(exprs(e))
  e <- e[keep,]

  design <- model.matrix( ~condition, coldt)
  dgel <- DGEList(exprs(e))
  dgel <- calcNormFactors(dgel)
  # current recommendation (Sep 2016) according to vignette:
  dgel <- estimateDisp(dgel, design)
  edger.fit <- glmFit(dgel, design)
  edger.lrt <- glmLRT(edger.fit)
  pvals <- edger.lrt$table$PValue

  padj[keep] <- p.adjust(pvals,method="BH")
  padj[is.na(padj)] <- 1

  pval[keep] <- edger.lrt$table$PValue
  pval[is.na(pval)] <- 1
  logFc[keep] <- edger.lrt$table$logFC

  # combine results
  res <- data.frame(padj = padj, pvalue = pval, log2FoldChange = logFc, row.names = gene.names)
  res
}

# Genewise Negative Binomial Generalized Linear Models With Quasi-Likelihood Tests
runEdgeRQL <- function(e, coldt) {
  padj <- rep(NA, nrow(e))
  pval <- rep(NA, nrow(e))
  logFc <- rep(NA, nrow(e))

  keep <- cpmFilter(e, cpm=TRUE)
  gene.names <- rownames(exprs(e))
  e <- e[keep,]

  design <- model.matrix( ~condition, coldt)
  dgel <- DGEList(exprs(e))
  dgel <- calcNormFactors(dgel)
  # current recommendation (Sep 2016) according to vignette:
  dgel <- estimateDisp(dgel, design)
  edger.fit <- glmQLFit(dgel,design)
  edger.qlf <- glmQLFTest(edger.fit)
  pvals <- edger.qlf$table$PValue

  padj[keep] <- p.adjust(pvals,method="BH")
  padj[is.na(padj)] <- 1

  pval[keep] <- edger.qlf$table$PValue
  pval[is.na(pval)] <- 1
  logFc[keep] <- edger.qlf$table$logFC

  # combine results
  res <- data.frame(ql.padj = padj, ql.pvalue = pval, log2FoldChange = logFc, row.names = gene.names)
  res
}


contrast <- snakemake@params[["contrast"]]
cts <- read.table(snakemake@input[["counts"]], header=TRUE, row.names="gene")
coldata <- read.table(snakemake@input[["samples"]], header=TRUE, row.names="sample")
colnames(cts) <- rownames(coldata)

# use only samples for the requested contrast
# run edgeR
idx <- as.character(coldata[, "condition"]) %in% contrast
edata   <-new("ExpressionSet", exprs=as.matrix(cts[, idx]))
cdt <- data.frame(condition = coldata[idx,'condition'], row.names = rownames(coldata)[idx])
cdt <- droplevels(cdt) # make sure the design matrix contains only used levels

res <- runEdgeR(edata, cdt)
resQL <- runEdgeRQL(edata, cdt)
res <- cbind(res, resQL[,c('ql.padj','ql.pvalue')])

write.table(as.data.frame(res), file=snakemake@output[["table"]])
