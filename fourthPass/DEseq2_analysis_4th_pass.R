##======================================================================
##  File  : DEseq2_analysis_4th_pass.R
##  Author: Jianying Li
##  Credit: https://www.huber.embl.de/users/klaus/Teaching/DESeq2Predoc2014.html#inspection-and-correction-of-pvalues
##  Credit: https://lashlock.github.io/compbio/R_presentation.html 
##  Comment: focus on mouse study
##======================================================================
source("x:/project2019/RNAseqProj/scripts/R_help_scripts/RNAseqHelpers.R")
setwd("x:/project2019/RNAseqProj/results/fourthPassAnalysis/")



##=================================
##  Consistent with third pass
##=================================

load("X:/project2019/RNAseqProj/results/Salmon/mouse_STAR_matrix_cnt.rda")
d.cleaned = mouseSTAR$matrix
metaData  = mouseSTAR$meta

d.cleaned.full = d.cleaned 
metaData.full = metaData

d.cleaned = d.cleaned.full[,-c(3,6)]
metaData = metaData.full [-c(2,5), ]


##  Get two samples only for mouse study -- fourth pass

mouseSTAR = list("matrix" = d.cleaned, "meta" = metaData)
save (mouseSTAR , file = "X:/project2019/RNAseqProj/results/fourthPassAnalysis/mouseResults/mouse_STAR_matrix_cnt_no_sample2.rda")


##===============================================
##  Follow the new protocol and 
##  Excluding the sample 2 (different batch)
##===============================================

##  Step 1: Create an DESeq2 object
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
#10054  4395

res     = DESeq2Res
res.raw = DESeq2Res

DEGs.adjp01.fc1.raw <- subset(res.raw, padj<.01 & abs(log2FoldChange)>1)

DESeq2Res = res.raw

##  back to thow two plots!!

##=================================
##  Now, try to correct p-values
##=================================


hist(DESeq2Res$pvalue, col = "lavender", main = "Mouse SH vs AG - p.value distribution", xlab = "p-values")
dim(DESeq2Res)
#[1] 18259     6
DESeq2Res <- DESeq2Res[ !is.na(DESeq2Res$padj), ]
dim(DESeq2Res)
#[1] 14449     6
DESeq2Res <- DESeq2Res[ !is.na(DESeq2Res$pvalue), ]
dim(DESeq2Res)
#[1] 14449     6

plotMA(DESeq2Res)


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
#13973    476

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

with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="Mouse RNAseq Volcano plot -- TR_SH vs. TR_AG", xlim=c(-10,10)))
# Add colored points: blue if padj<0.01, red if log2FC>1 and padj<0.05)
with(subset(res, padj<.01 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(res, padj<.01 & abs(log2FoldChange)>1), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
mtext ("adjusted p-value < 0.01, log2FoldChange > 1")
dim(subset(res, padj<.01 & abs(log2FoldChange)>1))


DEGs.adjp01.fc1 <- subset(res, padj<.01 & abs(log2FoldChange)>1)

normalizedReadCountsAll = counts(DESeq2Table,normalized=TRUE)
dim(normalizedReadCountsAll)
row.names(normalizedReadCountsAll)
mouseSTARDESeq2 = list("model" = dds, "results" = res, "DEGs_adjp01_fc1_dat" = normalizedReadCountsAll[which(rownames(normalizedReadCountsAll) %in% rownames(DEGs.adjp01.fc1)),])

save (mouseSTARDESeq2 , file = "X:/project2019/RNAseqProj/results/fourthPassAnalysis/mouse_STAR_DESeq2_results_noSample2.rda")
write.table (normalizedReadCountsAll[which(rownames(normalizedReadCountsAll) %in% rownames(DEGs.adjp01.fc1)),] , file = "X:/project2019/RNAseqProj/results/fourthPassAnalysis/mouse_STAR_DESeq2_DEGs_noSample2.txt", sep = "\t", row.names = TRUE, col.names = NA )



length(which(rownames(DEGs.adjp01.fc1.raw)  %in% rownames(DEGs.adjp01.fc1)))

data <-   normalizedReadCountsAll[which(rownames(normalizedReadCountsAll) %in% rownames(DEGs.adjp01.fc1)),]     
data <-   normalizedReadCountsAll[which(rownames(normalizedReadCountsAll) %in% rownames(DEGs.adjp01.fc1.raw)),] 
heatmap.2(as.matrix(data),   scale="row", key=T, keysize=1.5,
          density.info="none", trace="none",cexCol=0.9, labRow=NA)


dim(data)

heatmap.2(as.matrix(d.cleaned[-1]),   scale="row", key=T, keysize=1.5,
          density.info="none", trace="none",cexCol=0.9, labRow=NA)

##  There could be potentially batch effect with SH_2 or AG_2

##===============================================

