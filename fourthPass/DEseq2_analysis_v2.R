##======================================================================
##  File    : DEseq2_analysis_v2.R
##  Comment : with focus on human samples only
##  History : DEseq2_analysis.R
##  Author  : Jianying Li
##  Credit  : https://lashlock.github.io/compbio/R_presentation.html 
##======================================================================
source("x:/R-project/customPackages/dataManipTools.R")
source("x://R-project/customPackages/arraySeqTools.R")

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


##========================================
##  Don't nee it if don't work with it
##========================================

#countData <- read.csv("X:/project2019/RNAseqProj/results/counts/hg38_RNAseq_STAR_FeatureCounts_matrix.csv", header = TRUE)
#d.cleaned <- cleanMatrixData (countData)
#colnames(d.cleaned)
#dim(d.cleaned)
# [1] "EntrezGene" "K70_1"      "K70_2"      "K70_3"      "K70MD_1"    "K70MD_2"    "K70MD_3"    "K70P10_1"   "K70P10_2"   "K70P10_3"   "NKO70_1"    "NKO70_2"   
#[13] "NKO70_3" 

metaData <- read.csv("X:/project2019/RNAseqProj/meta/human_meta.csv", header = TRUE)
metaData 


load ("X:/project2019/RNAseqProj/results/firstPassAnalysis/K70ALL4_DEseq2.rda")
normalizedReadCountsAll = deseq2.results.K70.all4$norm.dat



##=======================================
##  Getting some heatmap -- clustering
##=======================================
library(gplots)

load ("X:/project2019/RNAseqProj/results/firstPassAnalysis/K70MD_DEseq2.rda")
res = deseq2.results.K70MD$results
K70MD.DEGs.high.stringency <- rownames(res[(which(res$padj <= 0.01& abs(res$log2FoldChange)>2)),])
K70MD.DEGs.low.stringency  <- rownames(res[(which(res$padj <= 0.01& abs(res$log2FoldChange)>1)),])


K70MD.DEGs.low.stringency  <- rownames(res[(which(res$padj <= 0.01& abs(res$log2FoldChange)>1)),])
length(K70MD.DEGs.low.stringency)
#[1] 6609
write.table(K70MD.DEGs , "X:/project2019/RNAseqProj/results/secondPassAnalysis/K70MD_DEGs.txt", sep="\t")


load("X:/project2019/RNAseqProj/results/firstPassAnalysis/P10_DEseq2.rda")
res = deseq2.results.K70P10$results
K70P10.DEGs.high.stringency  <- rownames(res[(which(res$padj <= 0.05& abs(res$log2FoldChange)>2)),])
K70P10.DEGs.low.stringency  <- rownames(res[(which(res$padj <= 0.05& abs(res$log2FoldChange)>1)),])

common.DEGs.K70MD.K70P10 <- intersect(K70MD.DEGs.low.stringency, K70P10.DEGs.low.stringency)
VennDiagram <- draw.two.list.venndigram(K70MD.DEGs.low.stringency, K70P10.DEGs.low.stringency, "K70MD" ,"P10")
common.DEGs.K70MD.K70P10.high<- intersect(K70MD.DEGs.high.stringency, K70P10.DEGs.high.stringency)

#write.table(common.DEGs.K70MD.K70P10.high, file = "x:/project2019/RNAseqProj/results/secondPassAnalysis/ESCC_196.txt", sep = "\t")

VennDiagram <- draw.two.list.venndigram(K70MD.DEGs.high.stringency, K70P10.DEGs.high.stringency, "K70MD" ,"P10")
VennDiagram <- draw.two.list.venndigram(K70MD.DEGs.high.stringency, K70P10.DEGs.low.stringency, "K70MD" ,"P10")
grid.draw(VennDiagram$figure)

VennDiagram.2 <- draw.three.list.venndigram(K70MD.DEGs.high.stringency,VennDiagram$union, K70P10.DEGs.high.stringency, "K70MD" ,"Common" ,"P10")
grid.draw(VennDiagram.2$figure)


load("X:/project2019/RNAseqProj/results/firstPassAnalysis/NKO70_DEseq2.rda")
res = deseq2.results.K70NKO70$results
K70NKO70.DEGs <- rownames(res[(which(res$padj <= 0.05& abs(res$log2FoldChange)>1)),])



length(common.DEGs.K70MD.K70P10.high)
#[1] 196
data <- normalizedReadCountsAll[which(rownames(normalizedReadCountsAll) %in% K70MD.DEGs),]
data <- normalizedReadCountsAll[which(rownames(normalizedReadCountsAll) %in% K70P10.DEGs),]
data <- normalizedReadCountsAll[which(rownames(normalizedReadCountsAll) %in% K70NKO70.DEGs),]
data <- normalizedReadCountsAll[which(rownames(normalizedReadCountsAll) %in% common.DEGs.K70MD.K70P10),]
data <- normalizedReadCountsAll[which(rownames(normalizedReadCountsAll) %in% common.DEGs.K70MD.K70P10.high),]



heatmap.2(as.matrix(data),   scale="row", key=T, keysize=1.5,
          density.info="none", trace="none",cexCol=0.9, labRow=NA)

write.table (data, file = "X:/project2019/RNAseqProj/results/secondPassAnalysis/dat_common_DEGs_K70MD_K70P10.txt",row.names = TRUE,col.names = NA, sep="\t")
write.table (data, file = "X:/project2019/RNAseqProj/results/secondPassAnalysis/dat_common_DEGs_K70MD_K70P10_high.txt",row.names = TRUE,col.names = NA, sep="\t")

##=============================================
##  Getting a venn-diagram of three gene sets
##=============================================
source("x:/R-project/customPackages/plotTools.R")
source("x:/R-project/customPackages/dataManipTools.R")
VennDiagram <- draw.three.list.venndigram(K70MD.DEGs, K70P10.DEGs, K70NKO70.DEGs , "K70MD" ,"P10", "NKO70")
grid.draw(VennDiagram$figure)



intersect.DEGs <- VennDiagram$union
data <- normalizedReadCountsAll[which(rownames(normalizedReadCountsAll) %in% intersect.DEGs ),]


heatmap.2(as.matrix(data),   scale="row", key=T, keysize=1.5,
          density.info="none", trace="none",cexCol=0.9, labRow=NA)

##  second pass analysis

##  first pass analysis

##===================================
##  K70 vs K70P10 (human)
##===================================
dds <- DESeqDataSetFromMatrix(countData=cleanMatrixData(d.cleaned [,c(1:4,8:10)]), 
                              colData=metaData [c(1:3,7:9),] , 
                              design=~treatment, tidy = TRUE)
dds <- DESeq(dds)
res <- results(dds)
dim(res[(which(res$padj <= 0.05 )),])

write.csv(res[(which(res$padj <= 0.05 )),], file = "x:/project2019/RNAseqProj/results/DE/human_K70P10_DGE_adjp05.csv")

png ("x:/project2019/RNAseqProj/results/DE/human_K70P10_DGE_adjp05.png")
par(mfrow=c(1,1))
# Make a basic volcano plot
with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="Human RNAseq Volcano plot", xlim=c(-3,3)))
# Add colored points: blue if padj<0.01, red if log2FC>1 and padj<0.05)
with(subset(res, padj<.05 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(res, padj<.05 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
mtext ("K70 vs K70P10")
dev.off()


##  first pass analysis

with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="Human RNAseq Volcano plot -- K70 vs K70P10", xlim=c(-3,3)))
# Add colored points: blue if padj<0.01, red if log2FC>1 and padj<0.05)
with(subset(res, padj<.05 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(res, padj<.05 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
mtext ("p-value < 0.05, foldChance > 2")
dim(subset(res, padj<.05 & abs(log2FoldChange)>2))

with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="Human RNAseq Volcano plot -- K70 vs K70P10", xlim=c(-3,3)))
# Add colored points: blue if padj<0.01, red if log2FC>1 and padj<0.05)
with(subset(res, padj<.01 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(res, padj<.01 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
mtext ("p-value < 0.01, foldChance > 2")

dim(subset(res, padj<.01 & abs(log2FoldChange)>2))
##  first pass analysis


##===================================
##  K70 vs NKO70 (human)
##===================================
dds <- DESeqDataSetFromMatrix(countData=cleanMatrixData(d.cleaned [,c(1:4,11:13)]), 
                              colData=metaData [c(1:3,10:12),] , 
                              design=~treatment, tidy = TRUE)
dds <- DESeq(dds)
res <- results(dds)
dim(res[(which(res$padj <= 0.05 )),])

write.csv(res[(which(res$padj <= 0.05 )),], file = "x:/project2019/RNAseqProj/results/DE/human_NKO70_DGE_adjp05.csv")


png ("x:/project2019/RNAseqProj/results/DE/human_NKO70_DGE_adjp05.png")
par(mfrow=c(1,1))
# Make a basic volcano plot
with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="Human RNAseq Volcano plot", xlim=c(-3,3)))
# Add colored points: blue if padj<0.01, red if log2FC>1 and padj<0.05)
with(subset(res, padj<.05 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(res, padj<.05 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
mtext ("K70 vs NKO70")
dev.off()


##  first pass analysis

with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="Human RNAseq Volcano plot -- K70 vs NKO70", xlim=c(-3,3)))
# Add colored points: blue if padj<0.01, red if log2FC>1 and padj<0.05)
with(subset(res, padj<.05 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(res, padj<.05 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
mtext ("p-value < 0.05, foldChance > 2")
dim(subset(res, padj<.05 & abs(log2FoldChange)>2))
    
with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="Human RNAseq Volcano plot -- K70 vs NKO70", xlim=c(-3,3)))
# Add colored points: blue if padj<0.01, red if log2FC>1 and padj<0.05)
with(subset(res, padj<.01 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(res, padj<.01 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
mtext ("p-value < 0.01, foldChance > 2")

dim(subset(res, padj<.01 & abs(log2FoldChange)>2))
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

write.table(normalizedReadCountsAll,file = "x:/project2019/RNAseqProj/results/DE/human_normalizedCnt.txt",row.names = TRUE,col.names = NA, sep="\t")


##===========================================

library(readxl)

PAM50 <- read_excel("x:/project2019/RNAseqProj/doc/PAM50.xlsx")
read.table(common.DEGs.K70MD.K70P10.high, file = "x:/project2019/RNAseqProj/results/secondPassAnalysis/ESCC_196.txt", sep = "\t")
pam50.gene.id <- PAM50$EntrezID

ESCC <- read_excel("x:/project2019/RNAseqProj/doc/Nrf2_genes_DHuang.xlsx")
ESCC_Nrf2 <- ESCC$`Gene Name`

temp.list <- c (common.DEGs.K70MD.K70P10.high, pam50.gene.id)
VennDiagram.3 <- draw.three.list.venndigram(common.DEGs.K70MD.K70P10.high, pam50.gene.id,ESCC_Nrf2, "ESCC196", "PAM50", "ESCCNrf2")
grid.draw(VennDiagram.3$figure)
