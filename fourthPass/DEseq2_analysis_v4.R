##======================================================================
##  File  : DEseq2_analysis_v4.R
##  Author: Jianying Li
##  Credit: https://www.huber.embl.de/users/klaus/Teaching/DESeq2Predoc2014.html#inspection-and-correction-of-pvalues
##  Credit: https://lashlock.github.io/compbio/R_presentation.html 
##  Comment: focus on mouse study
##======================================================================
source("x:/project2019/RNAseqProj/scripts/R_help_scripts/RNAseqHelpers.R")
setwd("x:/project2019/RNAseqProj/results/thirdPassAnalysis/")


countData <- read.csv("X:/project2019/RNAseqProj/results/Salmon/mouse_Salmon_matrix.csv", header = TRUE)
d.cleaned <- cleanMatrixData (countData)
colnames(d.cleaned)
dim(d.cleaned)

dm = d.cleaned

dm [,c(2:7)] = ceiling(dm[,c(2:7)])
head(dm)


metaData <- read.csv("X:/project2019/RNAseqProj/meta/mouse_meta.csv", header = TRUE)
metaData 
mouseSalmon = list("matrix" = dm, "meta" = metaData)

save (mouseSalmon, file = "X:/project2019/RNAseqProj/results/Salmon/mouse_Salmon_matrix_discrete.rda")

mouseSalmon$matrix[which(mouseSalmon$matrix$Ensemble_ID %in% "ENSMUST00000068106.4"),]
mouseSalmon$matrix[which(mouseSalmon$matrix$Ensemble_ID %in% "ENSMUST00000055562.2"),]


##  TR_AG vs TR_SH (mouse)
dds <- DESeqDataSetFromMatrix(countData=dm , 
                              colData=metaData , 
                              design=~treatment, tidy = TRUE)
dds <- DESeq(dds)
res <- results(dds)

mouseSalmonDESeq2 = list("model" = dds, "results" = res)
save (mouseSalmonDESeq2, file = "X:/project2019/RNAseqProj/results/Salmon/mouse_Salmon_DESeq2.rda")


mouseSalmonDESeq2$results[which(row.names(mouseSalmonDESeq2$results) %in% "ENSMUST00000068106.4"),]
mouseSalmonDESeq2$results[which(row.names(mouseSalmonDESeq2$results) %in% "ENSMUST00000055562.2"),]
##=========================================
##  DESeq2 is NOT appropriate for Salmon
##=========================================



##======================
## STAR/FeatureCounts
##======================

countData <- read.csv("X:/project2019/RNAseqProj/results/counts/mm10_RNAseq_STAR_FeatureCounts_matrix.csv", header = TRUE)
d.cleaned <- cleanMatrixData (countData)
dim(d.cleaned)
colnames(d.cleaned)

metaData <- read.csv("X:/project2019/RNAseqProj/meta/mouse_meta.csv", header = TRUE)
metaData 

mouseSTAR = list("matrix" = d.cleaned, "meta" = metaData)

save (mouseSTAR , file = "X:/project2019/RNAseqProj/results/Salmon/mouse_STAR_matrix_cnt.rda")



##  TR_AG vs TR_SH (mouse)
dds <- DESeqDataSetFromMatrix(countData=d.cleaned , 
                              colData=metaData , 
                              design=~treatment, tidy = TRUE)
dds <- DESeq(dds)
res <- results(dds)


##=================================
##  Consistent with third pass
##=================================

load("X:/project2019/RNAseqProj/results/Salmon/mouse_STAR_matrix_cnt.rda")
d.cleaned = mouseSTAR$matrix
metaData  = mouseSTAR$meta

##  TR_AG vs TR_SH (mouse)
DESeq2Table <- DESeqDataSetFromMatrix(countData=d.cleaned , 
                              colData=metaData , 
                              design=~treatment, tidy = TRUE)

##===================================
##  Getting the DEGs step-by-step
##  Turn out it is exactly same!!
##===================================
DESeq2Table <- estimateSizeFactors(DESeq2Table)
sizeFactors(DESeq2Table)

#  head(assay(DESeq2Table))
#  colSums(assay(DESeq2Table))
#  colData(DESeq2Table)
#  rowData(DESeq2Table)
#  mcols(rowData(DESeq2Table))

colData(DESeq2Table)$treatment 
DESeq2Table <- estimateDispersions(DESeq2Table)
plotDispEsts(DESeq2Table)
DESeq2Table <-  nbinomWaldTest(DESeq2Table)
DESeq2Res <- results(DESeq2Table, pAdjustMethod = "BH")
table(DESeq2Res$padj < 0.01)
#FALSE  TRUE 
#10506 10570


res = DESeq2Res

##  back to thow two plots!!

##=================================
##  Now, try to correct p-values
##=================================


hist(DESeq2Res$pvalue, col = "lavender", main = "K70MD vs K7", xlab = "p-values")
dim(DESeq2Res)
#[1] 18259     6
DESeq2Res <- DESeq2Res[ !is.na(DESeq2Res$padj), ]
dim(DESeq2Res)
#[1] 13453     6
DESeq2Res <- DESeq2Res[ !is.na(DESeq2Res$pvalue), ]
dim(DESeq2Res)
#[1] 13453     6


##==================================
##  try to adjust pvalues
##==================================

DESeq2Res <- DESeq2Res[, -which(names(DESeq2Res) == "padj")]
FDR.DESeq2Res <- fdrtool(DESeq2Res$stat, statistic= "normal", plot = T)
FDR.DESeq2Res$param[1, "sd"]
hist(FDR.DESeq2Res$pval, col = "royalblue4", 
     main = "K70MD vs K70 correct null model", xlab = "CORRECTED p-values")

DESeq2Res[,"padj"]  <- p.adjust(FDR.DESeq2Res$pval, method = "BH")


table(DESeq2Res[,"padj"] < 0.01)
#FALSE  TRUE 
#13370    83 

res = DESeq2Res

##=================================
##  Consistent with third pass
##  END
##=================================




d.cleaned[which(d.cleaned$Geneid == "Hoxc12"),]
d.cleaned[which(d.cleaned$Geneid == "6330419J24Rik"),]


##  first pass analysis



with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="Mouse RNAseq Volcano plot -- TR_SH vs. TR_AG", xlim=c(-10,10)))

# Add colored points: blue if padj<0.01, red if log2FC>1 and padj<0.05)
with(subset(res, padj<.05 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(res, padj<.05 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
mtext ("adjusted p-value < 0.05, log2FoldChange > 2")
dim(subset(res, padj<.05 & abs(log2FoldChange)>2))

with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="Mouse RNAseq Volcano plot -- TR_SH vs. TR_AG", xlim=c(-10,10)))
# Add colored points: blue if padj<0.01, red if log2FC>1 and padj<0.05)
with(subset(res, padj<.01 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(res, padj<.01 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
mtext ("adjusted p-value < 0.01, log2FoldChange > 2")
dim(subset(res, padj<.01 & abs(log2FoldChange)>2))

DEGs.adjp01.fc1 <- subset(res, padj<.01 & abs(log2FoldChange)>1)

normalizedReadCountsAll = counts(DESeq2Table,normalized=TRUE)
dim(normalizedReadCountsAll)
row.names(normalizedReadCountsAll)
mouseSTARDESeq2 = list("model" = dds, "results" = res, "DEGs_adjp01_fc1_dat" = normalizedReadCountsAll[which(rownames(normalizedReadCountsAll) %in% rownames(DEGs.adjp01.fc1)),])

save (mouseSTARDESeq2 , file = "X:/project2019/RNAseqProj/results/fourthPassAnalysis/mouse_STAR_DESeq2_results.rda")
write.table (normalizedReadCountsAll[which(rownames(normalizedReadCountsAll) %in% rownames(DEGs.adjp01.fc1)),] , file = "X:/project2019/RNAseqProj/results/fourthPassAnalysis/mouse_STAR_DESeq2_DEGs.txt", sep = "\t", row.names = TRUE, col.names = NA )

data <-   normalizedReadCountsAll[which(rownames(normalizedReadCountsAll) %in% rownames(DEGs.adjp01.fc1)),]                
heatmap.2(as.matrix(data),   scale="row", key=T, keysize=1.5,
          density.info="none", trace="none",cexCol=0.9, labRow=NA)


dim(data)

heatmap.2(as.matrix(d.cleaned[-1]),   scale="row", key=T, keysize=1.5,
          density.info="none", trace="none",cexCol=0.9, labRow=NA)

##  There could be potentially batch effect with SH_2 or AG_2

##===============================================