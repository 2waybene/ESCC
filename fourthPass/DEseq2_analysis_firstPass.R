##======================================================================
##  File  : DEseq2_analysis_firstPass.R
##  Author: Jianying Li
##  Credit: https://lashlock.github.io/compbio/R_presentation.html 
##======================================================================


cleanMatrixData <- function (dmIN, idColumn = 1)
{
  require (plyr)
  ##  clean data with 0 cross all samples
  t <- dmIN[-(which(rowSums(dmIN[,-idColumn]) == 0 )),]
  dim(t)
  
  ##  Filtering out row with ID as "?"
  if (length(which(t[,idColumn] == "?")) != 0 )
  {
    t.2 <- t[-(which(t[,idColumn] == "?")),]
    dim(t.2)
    t.3 <- ddply(t.2,colnames(t.2)[idColumn],numcolwise(sum))
    dim(t.3)
  }else{
    t.3 <- ddply(t,colnames(t)[idColumn],numcolwise(sum))
  }
  ##   remove duplicate entries by taking the sum
  return(t.3)
}
##===================================================================

library( "DESeq2" )
library(ggplot2)


countData <- read.csv("X:/project2019/RNAseqProj/results/counts/hg38_RNAseq_STAR_FeatureCounts_matrix.csv", header = TRUE)
d.cleaned <- cleanMatrixData (countData)
colnames(d.cleaned)
dim(d.cleaned)
# [1] "EntrezGene" "K70_1"      "K70_2"      "K70_3"      "K70MD_1"    "K70MD_2"    "K70MD_3"    "K70P10_1"   "K70P10_2"   "K70P10_3"   "NKO70_1"    "NKO70_2"   
#[13] "NKO70_3" 

metaData <- read.csv("X:/project2019/RNAseqProj/meta/human_meta.csv", header = TRUE)
metaData 

##==========================
##  K70 vs K70MD (human)
##=========================
dds <- DESeqDataSetFromMatrix(countData=cleanMatrixData(d.cleaned [,c(1:7)]), 
                              colData=metaData [c(1:6),] , 
                              design=~treatment, tidy = TRUE)
dds <- DESeq(dds)
res <- results(dds)
dim(res[(which(res$padj <= 0.05 )),])

#write.csv(res[(which(res$padj <= 0.05 )),], file = "x:/project2019/RNAseqProj/results/DE/human_K70MD_DGE_adjp05.csv")


##  first pass analysis

with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="Human RNAseq Volcano plot -- K70 vs K70MD", xlim=c(-3,3)))
# Add colored points: blue if padj<0.01, red if log2FC>1 and padj<0.05)
with(subset(res, padj<.05 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(res, padj<.05 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
mtext ("adjusted p-value < 0.05, log2FoldChange > 2")
dim(subset(res, padj<.05 & abs(log2FoldChange)>2))


with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="Human RNAseq Volcano plot -- K70 vs K70MD", xlim=c(-3,3)))
# Add colored points: blue if padj<0.01, red if log2FC>1 and padj<0.05)
with(subset(res, padj<.01 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(res, padj<.01 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
mtext ("adjusted p-value < 0.01, log2FoldChange > 2")
dim(subset(res, padj<.01 & abs(log2FoldChange)>2))
##  first pass analysis

deseq2.results.K70MD <- list (dat = cleanMatrixData(d.cleaned [,c(1:7)]), model = dds, results = res)
save (deseq2.results.K70MD, file = "X:/project2019/RNAseqProj/results/firstPassAnalysis/K70MD_DEseq2.rda")
##===================================
##  K70 vs K70P10 (human)
##===================================
dds <- DESeqDataSetFromMatrix(countData=cleanMatrixData(d.cleaned [,c(1:4,8:10)]), 
                              colData=metaData [c(1:3,7:9),] , 
                              design=~treatment, tidy = TRUE)
dds <- DESeq(dds)
res <- results(dds)
dim(res[(which(res$padj <= 0.05 )),])




##  first pass analysis

with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="Human RNAseq Volcano plot -- K70 vs K70P10", xlim=c(-3,3)))
# Add colored points: blue if padj<0.01, red if log2FC>1 and padj<0.05)
with(subset(res, padj<.05 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(res, padj<.05 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
mtext ("adjusted p-value < 0.05, log2FoldChange > 2")
dim(subset(res, padj<.05 & abs(log2FoldChange)>2))

with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="Human RNAseq Volcano plot -- K70 vs K70P10", xlim=c(-3,3)))
# Add colored points: blue if padj<0.01, red if log2FC>1 and padj<0.05)
with(subset(res, padj<.01 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(res, padj<.01 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
mtext ("adjusted p-value < 0.01, log2FoldChange > 2")

dim(subset(res, padj<.01 & abs(log2FoldChange)>2))
##  first pass analysis

deseq2.results.K70P10 <- list (dat = cleanMatrixData(d.cleaned [,c(1:7)]), model = dds, results = res)
save (deseq2.results.K70P10, file = "X:/project2019/RNAseqProj/results/firstPassAnalysis/P10_DEseq2.rda")
##===================================
##  K70 vs NKO70 (human)
##===================================
dds <- DESeqDataSetFromMatrix(countData=cleanMatrixData(d.cleaned [,c(1:4,11:13)]), 
                              colData=metaData [c(1:3,10:12),] , 
                              design=~treatment, tidy = TRUE)
dds <- DESeq(dds)
res <- results(dds)
dim(res[(which(res$padj <= 0.05 )),])



##  first pass analysis

with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="Human RNAseq Volcano plot -- K70 vs NKO70", xlim=c(-3,3)))
# Add colored points: blue if padj<0.01, red if log2FC>1 and padj<0.05)
with(subset(res, padj<.05 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(res, padj<.05 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
mtext ("adjusted p-value < 0.05, log2FoldChange > 2")
dim(subset(res, padj<.05 & abs(log2FoldChange)>2))



with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="Human RNAseq Volcano plot -- K70 vs NKO70", xlim=c(-3,3)))
# Add colored points: blue if padj<0.01, red if log2FC>1 and padj<0.05)
with(subset(res, padj<.05 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(res, padj<.05 & abs(log2FoldChange)>1), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
mtext ("adjusted p-value < 0.05, log2FoldChange > 1")
dim(subset(res, padj<.05 & abs(log2FoldChange)>1))

    
with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="Human RNAseq Volcano plot -- K70 vs NKO70", xlim=c(-3,3)))
# Add colored points: blue if padj<0.01, red if log2FC>1 and padj<0.05)
with(subset(res, padj<.01 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(res, padj<.01 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
mtext ("adjusted p-value < 0.01, log2FoldChange > 2")
dim(subset(res, padj<.01 & abs(log2FoldChange)>2))

with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="Human RNAseq Volcano plot -- K70 vs NKO70", xlim=c(-3,3)))
# Add colored points: blue if padj<0.01, red if log2FC>1 and padj<0.05)
with(subset(res, padj<.01 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(res, padj<.01 & abs(log2FoldChange)>1), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
mtext ("adjusted p-value < 0.01, log2FoldChange > 1")
dim(subset(res, padj<.01 & abs(log2FoldChange)>1))


deseq2.results.K70NKO70 <- list (dat = cleanMatrixData(d.cleaned [,c(1:7)]), model = dds, results = res)
save (deseq2.results.K70NKO70, file = "X:/project2019/RNAseqProj/results/firstPassAnalysis/NKO70_DEseq2.rda")
##  first pass analysis


##===================================
##  all four groups (human)
##===================================

dds <- DESeqDataSetFromMatrix(countData=d.cleaned , 
                              colData=metaData  , 
                              design=~treatment, tidy = TRUE)
dds <- DESeq(dds)
res <- results(dds)
dim(res[(which(res$padj <= 0.05 )),])


write.csv(res[(which(res$padj <= 0.05 )),], file = "x:/project2019/RNAseqProj/results/DE/human_4_groups_DGE_adjp05.csv")

png ("x:/project2019/RNAseqProj/results/DE/human_4_groups_DGE_adjp05.png")
par(mfrow=c(1,1))
# Make a basic volcano plot
with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="Human RNAseq Volcano plot", xlim=c(-3,3)))
# Add colored points: blue if padj<0.01, red if log2FC>1 and padj<0.05)
with(subset(res, padj<.05 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(res, padj<.05 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
mtext ("All four groups")
dev.off()

normalizedReadCountsAll = counts(dds,normalized=TRUE)
dim(normalizedReadCountsAll)
head(normalizedReadCountsAll)

write.table(normalizedReadCountsAll,file = "x:/project2019/RNAseqProj/results/firstPassAnalysis/human_normalizedCnt.txt",row.names = TRUE,col.names = NA, sep="\t")

deseq2.results.K70.all4 <- list (dat = d.cleaned, model = dds, results = res, norm.dat = normalizedReadCountsAll)
save (deseq2.results.K70.all4, file = "X:/project2019/RNAseqProj/results/firstPassAnalysis/K70ALL4_DEseq2.rda")

##==========================
##  Mouse study
##==========================

countData <- read.csv("X:/project2019/RNAseqProj/results/counts/mm10_RNAseq_STAR_FeatureCounts_matrix.csv", header = TRUE)
d.cleaned <- cleanMatrixData (countData)
dim(d.cleaned)
colnames(d.cleaned)

metaData <- read.csv("X:/project2019/RNAseqProj/meta/mouse_meta.csv", header = TRUE)
metaData 

##  TR_AG vs TR_SH (mouse)
dds <- DESeqDataSetFromMatrix(countData=d.cleaned , 
                              colData=metaData , 
                              design=~treatment, tidy = TRUE)
dds <- DESeq(dds)
res <- results(dds)


##  first pass analysis

with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="Mouse RNAseq Volcano plot -- TR_SH vs. TR_AG", xlim=c(-3,3)))
# Add colored points: blue if padj<0.01, red if log2FC>1 and padj<0.05)
with(subset(res, padj<.05 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(res, padj<.05 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
mtext ("adjusted p-value < 0.05, log2FoldChange > 2")
dim(subset(res, padj<.05 & abs(log2FoldChange)>2))

with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="Mouse RNAseq Volcano plot -- TR_SH vs. TR_AG", xlim=c(-3,3)))
# Add colored points: blue if padj<0.01, red if log2FC>1 and padj<0.05)
with(subset(res, padj<.01 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(res, padj<.01 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
mtext ("adjusted p-value < 0.01, log2FoldChange > 2")
dim(subset(res, padj<.01 & abs(log2FoldChange)>2))


with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="Mouse RNAseq Volcano plot -- TR_SH vs. TR_AG", xlim=c(-3,3)))
# Add colored points: blue if padj<0.01, red if log2FC>1 and padj<0.05)
with(subset(res, padj<.05 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(res, padj<.05 & abs(log2FoldChange)>1), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
mtext ("adjusted p-value < 0.05, log2FoldChange > 1")
dim(subset(res, padj<.05 & abs(log2FoldChange)>1))

with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="Mouse RNAseq Volcano plot -- TR_SH vs. TR_AG", xlim=c(-3,3)))
# Add colored points: blue if padj<0.01, red if log2FC>1 and padj<0.05)
with(subset(res, padj<.01 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(res, padj<.01 & abs(log2FoldChange)>1), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
mtext ("adjusted p-value < 0.01, log2FoldChange > 1")
dim(subset(res, padj<.01 & abs(log2FoldChange)>1))

##  first pass analysis



