library(ggplot2)
library(DESeq2)
library(apeglm)
library(ggrepel)
library(dplyr)
library(org.Mm.eg.db)
library(PCAtools)
library(GEOquery)
library(fgsea)


data <- read.table("GSE186169_COUNTS_genes_LUNGCANCER_16.txt", sep="\t", header = 1)
head(data)



mapping <- data[, c("id_gene", "gene_name")]
rownames(mapping) <- mapping[, 1, drop=T]
head(mapping)



counts <- data[, 5:ncol(data)]
rownames(counts) <- mapping[, 1, drop=T]
head(counts)



annotations <- getGEO("GSE186169")[[1]]
head(pData(annotations))


pdata <- pData(annotations)[, c("title",
                                "description",
                                "age:ch1",
                                "cell line:ch1",
                                "Sex:ch1",
                                "tissue:ch1",
                                "treatment:ch1")]


rownames(pdata) <- pdata$description
colnames(pdata) <- c("Title", "Description", "Age", "Cells", "Sex", "Tissue", "Treatment")

pdata$Treatment <- as.character(pdata$Treatment)
pdata[is.na(pdata)] <- "No"
pdata$Treatment <- factor(pdata$Treatment, levels=c("No", "IFNg"))


dds <- DESeqDataSetFromMatrix(countData = counts[rowSums(counts) > 10, ],
                              colData = pdata,
                              design= ~ Treatment + Cells)
dds <- DESeq(dds)
resultsNames(dds)


vst <- varianceStabilizingTransformation(dds)
plotPCA(vst, intgroup=c("Cells"), ntop=nrow(vst)) + theme_bw() +
  theme(aspect.ratio = 1)

plotPCA(vst, intgroup=c("Treatment"), ntop=nrow(vst)) + theme_bw() +
  theme(aspect.ratio = 1)


res <- lfcShrink(dds, coef="Treatment_IFNg_vs_No", type="apeglm")
res$Gene.symbol <- mapping[rownames(res), "gene_name"]
head(res)


resDF <- as.data.frame(res)
ggplot(resDF, aes(x=log2FoldChange, y=-log10(padj), color=padj < 0.05)) +
  geom_point() + theme_bw() + scale_color_manual(values=c("black", "red")) +
  geom_text_repel(data=resDF %>% dplyr::filter(padj < 1e-20), aes(label=Gene.symbol), color="black") +
  xlim(c(-6, 6)) + xlab("Control     Log2FC    IFNg") + geom_vline(xintercept = 0, lty=2) +
  theme(aspect.ratio = 1)



upRegulatedGenes <- resDF %>% filter(padj < 0.01 & log2FoldChange > 0) %>% pull("Gene.symbol")
writeLines(upRegulatedGenes, "upregulated.txt")
head(upRegulatedGenes)


length(upRegulatedGenes)

dim(vst)


res <- results(dds, name="Treatment_IFNg_vs_No")
res$Gene.symbol <- mapping[rownames(res), "gene_name"]

stats <- res$stat
names(stats) <- res$Gene.symbol

load("keggSymbolHuman.rdata")
fgseaResults <- fgseaMultilevel(keggSymbolHuman, stats, minSize = 15, maxSize = 500)


ifn_genes <- readLines("ifn_gamma_genes.txt")
plotEnrichment(ifn_genes, stats)

topPathwaysUp <- fgseaResults[ES > 0, ][head(order(pval), n=5), pathway]
topPathwaysDown <- fgseaResults[ES < 0, ][head(order(pval), n=5), pathway]
topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
plotGseaTable(keggSymbolHuman[topPathways], stats, fgseaResults, gseaParam = 0.5,
              pathwayLabelStyle=list(size=10))

