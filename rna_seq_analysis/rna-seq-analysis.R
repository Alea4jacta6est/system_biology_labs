library(ggplot2)
library(DESeq2)
library(apeglm)
library(ggrepel)
library(dplyr)
library(org.Mm.eg.db)
library(PCAtools)


countFiles <- list.files("~/Code/sys_bio/Part_2/materials/GSE137888_RAW", full.names = T)
countFiles


counts <- lapply(countFiles, function(countsFile) {
  read.table(countsFile, sep="\t", header=1, skip = 1, row.names = 1, stringsAsFactors = F, comment.char = "")
})

head(counts[[1]])


counts <- lapply(counts, function(countsTable) countsTable[, "Count", drop=F])
counts <- do.call(cbind, counts)
colnames(counts) <- gsub(".*(SRR\\d+).*", "\\1", countFiles)
head(counts)


coldata <- data.frame(
  srr=gsub(".*(SRR\\d+).*", "\\1", countFiles),
  condition=gsub(".*(SRR\\d+)_(Control|Prdm16).*", "\\2", countFiles),
  row.names =gsub(".*(SRR\\d+).*", "\\1", countFiles)
)

write.table(coldata, file="~/Code/sys_bio/Part_2/materials/GSE137888_coldata.tsv", sep="\t", quote=F, col.names = NA)
coldata



dds <- DESeqDataSetFromMatrix(countData = counts[rowSums(counts) > 10, ],
                              colData = coldata,
                              design= ~ condition)
dds <- DESeq(dds)
resultsNames(dds) # lists the coefficients



rlog <- rlogTransformation(dds)
write.table(assay(rlog), file="~/Code/sys_bio/Part_2/materials/GSE137888_rlog.tsv", sep="\t", quote=F, col.names = NA)



vst <- varianceStabilizingTransformation(dds)
plotPCA(vst, intgroup=c("condition"), ntop=nrow(vst)) + theme_bw()



pcaData <- pca(assay(vst), metadata=coldata)
biplot(pcaData, colby="condition", legendPosition = "right")


plotloadings(pcaData) + theme_bw(base_size=8)



head(results(dds))


res <- lfcShrink(dds, coef="condition_Prdm16_vs_Control", type="apeglm")
head(res)


keytypes(org.Mm.eg.db)
res$Gene.symbol <- mapIds(org.Mm.eg.db, gsub("\\.\\d+", "", rownames(res)), column="SYMBOL", keytype="ENSEMBL")
head(res)


resDF <- as.data.frame(res)
volcanoPlot <- ggplot(resDF, aes(x=log2FoldChange, y=-log10(padj), color=padj < 0.05)) +
  geom_point() + theme_bw() + scale_color_manual(values=c("black", "red"))
volcanoPlot

volcanoPlot <- volcanoPlot +
  geom_text_repel(data=resDF %>% dplyr::filter(padj < 1e-20), aes(label=Gene.symbol), color="black")
volcanoPlot

volcanoPlot <- volcanoPlot +
  xlim(c(-12, 12)) +
  xlab("Prdm16     Log2FC    Control") +
  geom_vline(xintercept = 0, lty=2) +
  theme(aspect.ratio = 1)
volcanoPlot


library(fgsea)

deResults <- results(dds)
deResults$Gene.symbol <- mapIds(org.Mm.eg.db, gsub("\\.\\d+", "", rownames(deResults)), column="SYMBOL", "ENSEMBL")
stats <- deResults$stat
names(stats) <- deResults$Gene.symbol

load("~/Code/sys_bio/Part_2/materials/keggSymbolMouse.rdata")
fgseaResults <- fgseaMultilevel(keggSymbolMouse, stats, minSize = 15, maxSize = 500)
head(fgseaResults, 3)

topPathwaysUp <- fgseaResults[ES > 0, ][head(order(pval), n=8), pathway]
topPathwaysDown <- fgseaResults[ES < 0, ][head(order(pval), n=8), pathway]
topPathways <- c(topPathwaysUp, rev(topPathwaysDown))

plotGseaTable(keggSymbolMouse[topPathways], stats, fgseaResults, gseaParam = 0.5,
              pathwayLabelStyle=list(size=10))

