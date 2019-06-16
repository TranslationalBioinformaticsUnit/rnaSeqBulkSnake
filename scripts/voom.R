log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library("limma")
library("Biobase")
library("SummarizedExperiment")

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

runVoom <- function(e, coldt) {
  padj <- rep(NA, nrow(e))
  pval <- rep(NA, nrow(e))
  logFc <- rep(NA, nrow(e))

  gene.names <- rownames(exprs(e))
  keep <- cpmFilter(e, cpm=FALSE)
  e <- e[keep,]
  
  design <- model.matrix(~ condition, coldt)
  dgel <- DGEList(exprs(e))
  dgel <- calcNormFactors(dgel)
  v <- voom(dgel,design,plot=FALSE)
  fit <- lmFit(v,design)
  fit <- eBayes(fit)
  tt <- topTable(fit,coef=ncol(design),n=nrow(dgel),sort.by="none")
  pvals <- tt$P.Value 

  padj[keep] <- p.adjust(pvals,method="BH")
  padj[is.na(padj)] <- 1
  
  pval[keep] <- tt$P.Value
  pval[is.na(pval)] <- 1

  logFc[keep] <- tt$logFC

  # combine results
  res <- data.frame(padj = padj, pvalue = pval, log2FoldChange = logFc, row.names = gene.names)
  res
}

contrast <- snakemake@params[["contrast"]]
cts <- read.table(snakemake@input[["counts"]], header=TRUE, row.names="gene")
coldata <- read.table(snakemake@input[["samples"]], header=TRUE, row.names="sample")
colnames(cts) <- rownames(coldata)

# use only samples for the requested contrast
# run voom
idx <- as.character(coldata[, "condition"]) %in% contrast
edata   <-new("ExpressionSet", exprs=as.matrix(cts[, idx]))
cdt <- data.frame(condition = coldata[idx,'condition'], row.names = rownames(coldata)[idx])
cdt <- droplevels(cdt) # make sure the design matrix contains only used levels
res <- runVoom(edata, cdt)

write.table(as.data.frame(res), file=snakemake@output[["table"]])
