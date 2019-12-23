##======================================================================
##  File  : DEseq2_analysis_finalization.R
##  Author: Jianying Li
##  comment: to compare and finalize the DESeq method
##  Credit: https://www.huber.embl.de/users/klaus/Teaching/DESeq2Predoc2014.html#inspection-and-correction-of-pvalues
##  Credit: https://lashlock.github.io/compbio/R_presentation.html 
##======================================================================
source("x:/project2019/RNAseqProj/scripts/R_help_scripts/RNAseqHelpers.R")
setwd("x:/project2019/RNAseqProj/results/thirdPassAnalysis/")

##======================================
##  Based on the quick start, 
##  One can get DEGs easily
##  K70 vs K70MD (human)
##======================================

load ("X:/project2019/RNAseqProj/results/firstPassAnalysis/K70MD_DEseq2.rda")
res = deseq2.results.K70MD$results
table(res$padj < 0.01)

#FALSE  TRUE 
#10506 10570 

##  first pass analysis

with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="Human RNAseq Volcano plot -- K70 vs K70MD", xlim=c(-3,3)))
# Add colored points: blue if padj<0.01, red if log2FC>1 and padj<0.05)
with(subset(res, padj<.05 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(res, padj<.05 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
mtext ("adjusted p-value < 0.05, log2FoldChange > 2")
dim(subset(res, padj<.05 & abs(log2FoldChange)>2))
#[1] 2550    6

with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="Human RNAseq Volcano plot -- K70 vs K70MD", xlim=c(-3,3)))
# Add colored points: blue if padj<0.01, red if log2FC>1 and padj<0.05)
with(subset(res, padj<.01 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(res, padj<.01 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
mtext ("adjusted p-value < 0.01, log2FoldChange > 2")
dim(subset(res, padj<.01 & abs(log2FoldChange)>2))
#[1] 2017    6


##  Second pass analysis

K70MD.DEGs.high.stringency <- rownames(res[(which(res$padj <= 0.01& abs(res$log2FoldChange)>2)),])
K70MD.DEGs.low.stringency  <- rownames(res[(which(res$padj <= 0.01& abs(res$log2FoldChange)>1)),])

length(K70MD.DEGs.high.stringency)
# 2017
length(K70MD.DEGs.low.stringency)
#[1] 6609


##================================
##  testing Klaus' data
##  learing code
##  Third pass
##======================================
# load(url("http://www-huber.embl.de/users/klaus/geneCounts.RData"))

metaData <- read.csv("X:/project2019/RNAseqProj/meta/human_meta.csv", header = TRUE)
metaData 

##==========================
##  K70 vs K70MD (human)
##=========================
#countData <- read.csv("X:/project2019/RNAseqProj/results/counts/hg38_RNAseq_STAR_FeatureCounts_matrix.csv", header = TRUE)
#d.cleaned <- cleanMatrixData (countData)
#save (d.cleaned , file = "X:/project2019/RNAseqProj/results/counts/hg38_RNAseq_STAR_FeatureCounts_matrix_cleaned.rda")

load ("X:/project2019/RNAseqProj/results/counts/hg38_RNAseq_STAR_FeatureCounts_matrix_cleaned.rda")
dds <- DESeqDataSetFromMatrix(countData=cleanMatrixData(d.cleaned [,c(1:7)]), 
                              colData=metaData [c(1:6),] , 
                              design=~treatment, tidy = TRUE)

##===================================
##  Getting the DEGs step-by-step
##  Turn out it is exactly same!!
##===================================
DESeq2Table = dds 
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
#[1] 26790     6
DESeq2Res <- DESeq2Res[ !is.na(DESeq2Res$padj), ]
dim(DESeq2Res)
#[1] 21076     6
DESeq2Res <- DESeq2Res[ !is.na(DESeq2Res$pvalue), ]
dim(DESeq2Res)
#[1] 21076     6


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
#19306  1770 

res = DESeq2Res

with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="Human RNAseq Volcano plot -- K70 vs K70MD", xlim=c(-3,3)))
# Add colored points: blue if padj<0.01, red if log2FC>1 and padj<0.05)
with(subset(res, padj<.05 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(res, padj<.05 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
mtext ("adjusted p-value < 0.05, log2FoldChange > 2")
dim(subset(res, padj<.05 & abs(log2FoldChange)>2))
#[1] 579    6


with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="Human RNAseq Volcano plot -- K70 vs K70MD", xlim=c(-3,3)))
# Add colored points: blue if padj<0.01, red if log2FC>1 and padj<0.05)
with(subset(res, padj<.01 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(res, padj<.01 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
mtext ("adjusted p-value < 0.01, log2FoldChange > 2")
dim(subset(res, padj<.01 & abs(log2FoldChange)>2))
#[1] 447    6


with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="Human RNAseq Volcano plot -- K70 vs K70MD", xlim=c(-3,3)))
# Add colored points: blue if padj<0.01, red if log2FC>1 and padj<0.05)
with(subset(res, padj<.01 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(res, padj<.01 & abs(log2FoldChange)>1), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
mtext ("adjusted p-value < 0.01, log2FoldChange > 1")
dim(subset(res, padj<.01 & abs(log2FoldChange)>1))
#[1] 1727    6


#FALSE  TRUE 
#19306  1770 

plotMA(DESeq2Res)

sigGenes <- rownames(subset(DESeq2Res, padj < 0.01))
sigResults <- subset(DESeq2Res, padj < 0.01)
head(sigResults)
dim(sigResults)
length(sigResults$log2FoldChange > 1 | sigResults$log2FoldChange < -1) 


intersect (sigGenes , K70MD.DEGs.high.stringency)
length(intersect (sigGenes , K70MD.DEGs.low.stringency))
length(sigGenes)
length(K70MD.DEGs.low.stringency)


K70MD.DEGs.high.stringency <- rownames(res[(which(res$padj <= 0.01& abs(res$log2FoldChange)>2)),])
K70MD.DEGs.low.stringency  <- rownames(res[(which(res$padj <= 0.01& abs(res$log2FoldChange)>1)),])

res = DESeq2Res
##  first pass analysis

with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="Human RNAseq Volcano plot -- K70 vs K70MD", xlim=c(-3,3)))
# Add colored points: blue if padj<0.01, red if log2FC>1 and padj<0.05)
with(subset(res, padj<.05 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(res, padj<.05 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
mtext ("adjusted p-value < 0.05, log2FoldChange > 2")
dim(subset(res, padj<.05 & abs(log2FoldChange)>2))
rownames(subset(res, padj<.05 & abs(log2FoldChange)>2))

with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="Human RNAseq Volcano plot -- K70 vs K70MD", xlim=c(-3,3)))
# Add colored points: blue if padj<0.01, red if log2FC>1 and padj<0.05)
with(subset(res, padj<.01 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(res, padj<.01 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
mtext ("adjusted p-value < 0.01, log2FoldChange > 2")
dim(subset(res, padj<.01 & abs(log2FoldChange)>2))
##  first pass analysis


library(readxl)

PAM50 <- read_excel("x:/project2019/RNAseqProj/doc/PAM50.xlsx")
pam50.gene.id <- PAM50$EntrezID

ESCC <- read_excel("x:/project2019/RNAseqProj/doc/Nrf2_genes_DHuang.xlsx")
ESCC_Nrf2 <- ESCC$`Gene Name`


DEGs.K70MD <- rownames(subset(res, padj<.01 & abs(log2FoldChange)>2))

temp.list <- c (DEGs.K70MD, pam50.gene.id)
VennDiagram.3 <- draw.three.list.venndigram(DEGs.K70MD, pam50.gene.id,ESCC_Nrf2, "ESCC196", "PAM50", "ESCCNrf2")
grid.draw(VennDiagram.3$figure)

VennDiagram.3 <- draw.three.list.venndigram(DEGs.P10, pam50.gene.id,ESCC_Nrf2, "ESCC196", "PAM50", "ESCCNrf2")
VennDiagram.3 <- draw.two.list.venndigram(DEGs.P10, DEGs.K70MD, "P10", "MD")

#deseq2.results.K70MD <- list (dat = cleanMatrixData(d.cleaned [,c(1:7)]), model = dds, results = res)
#save (deseq2.results.K70MD, file = "X:/project2019/RNAseqProj/results/firstPassAnalysis/K70MD_DEseq2.rda")

##===================================
##  K70 vs K70P10 (human)
##===================================
##  exclude K70P10_3
#dds <- DESeqDataSetFromMatrix(countData=cleanMatrixData(d.cleaned [,c(1:4,8:10)]), 
# colData=metaData [c(1:3,7:9),] , 

dds <- DESeqDataSetFromMatrix(countData=cleanMatrixData(d.cleaned [,c(1:4,8:9)]), 
                              colData=metaData [c(1:3,7:8),] , 
                              design=~treatment, tidy = TRUE)
#dds <- DESeq(dds)
#res <- results(dds)
#dim(res[(which(res$padj <= 0.05 )),])



DESeq2Table = dds 
head(assay(DESeq2Table))
colSums(assay(DESeq2Table))
colData(DESeq2Table)
rowData(DESeq2Table)

mcols(rowData(DESeq2Table))


con <- as.character(colData(DESeq2Table)$treatment)
colData(DESeq2Table)$condition <- factor(con)

GeneCounts <- counts(DESeq2Table)
idx.nz <- apply(GeneCounts, 1, function(x) { all(x > 0)})
sum(idx.nz)

nz.counts <- subset(GeneCounts, idx.nz)
sam <- sample(dim(nz.counts)[1], 5)
nz.counts[sam, ]

colData(DESeq2Table)$treatment <- factor(ifelse(is.na(colData(DESeq2Table)$treatment),  "ctrl", "P10"), levels = c("ctrl", "P10"))

DESeq2Table <- estimateSizeFactors(DESeq2Table)
sizeFactors(DESeq2Table)

multidensity( counts(DESeq2Table, normalized = T)[idx.nz ,],
              xlab="mean counts", xlim=c(0, 1000))


multiecdf( counts(DESeq2Table, normalized = T)[idx.nz ,],
           xlab="mean counts", xlim=c(0, 1000))

pdf("pairwiseMAs_P10.pdf")
MA.idx = t(combn(1:5, 2))
for( i in  seq_along( MA.idx[,1])){ 
  MDPlot(counts(DESeq2Table, normalized = T)[idx.nz ,], 
         c(MA.idx[i,1],MA.idx[i,2]), 
         main = paste( colnames(DESeq2Table)[MA.idx[i,1]], " vs ",
                       colnames(DESeq2Table)[MA.idx[i,2]] ), ylim = c(-3,3))
}
dev.off()

rld <- rlogTransformation(DESeq2Table, blind=TRUE)


distsRL <- dist(t(assay(rld)))
mat <- as.matrix(distsRL)
#rownames(mat) <-  colData(rld)$treatment
#colnames(mat) <-  colData(rld)$sampleNO

hmcol <- colorRampPalette(brewer.pal(9, "Blues"))(255)
heatmap.2(mat, trace="none", col = rev(hmcol), margin=c(13, 13))



##    get a PCA analysis

ntop = 500
Pvars <- rowVars(assay(rld))
select <- order(Pvars, decreasing = TRUE)[seq_len(min(ntop, 
                                                      length(Pvars)))]
PCA <- prcomp(t(assay(rld)[select, ]), scale = F)
percentVar <- round(100*PCA$sdev^2/sum(PCA$sdev^2),1)


dataGG = data.frame(PC1 = PCA$x[,1], PC2 = PCA$x[,2], 
                    PC3 = PCA$x[,3], PC4 = PCA$x[,4], 
                    sampleNO =  colData(rld)$condition,
                    condition =  colData(rld)$condition)

(qplot(PC1, PC2, data = dataGG, color =  condition, 
       main = "PC1 vs PC2, top variable genes", size = I(6))
  + labs(x = paste0("PC1, VarExp:", round(percentVar[1],4)),
         y = paste0("PC2, VarExp:", round(percentVar[2],4)))
  + scale_colour_brewer(type="qual", palette=2)
)


##  No need for these
# outliers <- as.character(subset(colnames(DESeq2Table), dataGG$PC1 > 0))
# outliers

DESeq2Table <- estimateDispersions(DESeq2Table)
plotDispEsts(DESeq2Table)

DESeq2Table <-  nbinomWaldTest(DESeq2Table)
DESeq2Res <- results(DESeq2Table, pAdjustMethod = "BH")

table(DESeq2Res$padj < 0.01)

plot(metadata(DESeq2Res)$filterNumRej, type="b", xlab="quantiles of 'baseMean'",
     ylab="number of rejections")

hist(DESeq2Res$pvalue, col = "lavender", main = "K70 vs K70MD", xlab = "p-values")

DESeq2Res <- DESeq2Res[ !is.na(DESeq2Res$padj), ]
DESeq2Res <- DESeq2Res[ !is.na(DESeq2Res$pvalue), ]

DESeq2Res <- DESeq2Res[, -which(names(DESeq2Res) == "padj")]
FDR.DESeq2Res <- fdrtool(DESeq2Res$stat, statistic= "normal", plot = T)

FDR.DESeq2Res$param[1, "sd"]


DESeq2Res[,"padj"]  <- p.adjust(FDR.DESeq2Res$pval, method = "BH")
hist(FDR.DESeq2Res$pval, col = "royalblue4", 
     main = "K70 vs K70MD correct null model", xlab = "CORRECTED p-values")


table(DESeq2Res[,"padj"] < 0.01)

#FALSE  TRUE 
#19306  1770 

plotMA(DESeq2Res)

sigGenes <- rownames(subset(DESeq2Res, padj < 0.01))

sigResults <- subset(DESeq2Res, padj < 0.01)
head(sigResults)
dim(sigResults)

res = DESeq2Res
##  first pass analysis




##================================
##  Not working from here
##================================
anno <- AnnotationDbi::select(org.Mm.eg.db, 
                              keys=rownames(DESeq2Res), 
                              columns=c("SYMBOL", "GENENAME"),
                              keytype="SYMBOL")

anSig <- as.data.frame(subset(anno, SYMBOL %in% sigGenes))

sample_n(anSig, 5)


overallBaseMean <- as.matrix(DESeq2Res[, "baseMean", drop = F])

sig_idx <- match(anSig$ENSEMBL, rownames(overallBaseMean))

backG <- c()

for(i in sig_idx){
  ind <- genefinder(overallBaseMean, i, 10, method = "manhattan")[[1]]$indices
  backG <- c(backG, ind)
  
}

backG <- unique(backG)
backG <- rownames(overallBaseMean)[backG]

##================================================
##  Save the processed files for future use
##  Don't need to the following lines
##================================================
# countData <- read.csv("X:/project2019/RNAseqProj/results/counts/hg38_RNAseq_STAR_FeatureCounts_matrix.csv", header = TRUE)
# d.cleaned <- cleanMatrixData (countData)
# colnames(d.cleaned)
# dim(d.cleaned)
# [1] "EntrezGene" "K70_1"      "K70_2"      "K70_3"      "K70MD_1"    "K70MD_2"    "K70MD_3"    "K70P10_1"   "K70P10_2"   "K70P10_3"   "NKO70_1"    "NKO70_2"   
#[13] "NKO70_3" 

# metaData <- read.csv("X:/project2019/RNAseqProj/meta/human_meta.csv", header = TRUE)
# metaData 

# human_nrf2_dat <- list (dat = d.cleaned, meta = metaData)
# save (human_nrf2_dat, file  ="X:/project2019/RNAseqProj/results/counts/human_nrf2_dat.rda")
##================================================
##  Don't need to the following lines
##================================================

load("X:/project2019/RNAseqProj/results/counts/human_nrf2_dat.rda")
metaData = human_nrf2_dat$meta
d.cleaned = human_nrf2_dat$dat
metaData 
colnames(d.cleaned)
dim(d.cleaned)


##==========================
##  K70 vs K70MD (human)
##=========================
GeneCounts <- cleanMatrixData(d.cleaned [,c(1:7)])
head(GeneCounts)
idx.nz <- apply(GeneCounts, 1, function(x) { all(x > 0)})
sum(idx.nz)


##  create a DESeqDataSet object

multidensity( counts(DESeq2Table, normalized = T)[idx.nz ,],
              xlab="mean counts", xlim=c(0, 1000))

dds <- DESeqDataSetFromMatrix(countData=cleanMatrixData(d.cleaned [,c(1:7)]), 
                              colData=metaData [c(1:6),] , 
                              design=~treatment, tidy = TRUE)


counts(dds)
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

DEGs.K70MD <- subset(res, padj<.01 & abs(log2FoldChange)>2)
deseq2.results.K70MD <- list (dat = cleanMatrixData(d.cleaned [,c(1:7)]), model = dds, results = res)
save (deseq2.results.K70MD, file = "X:/project2019/RNAseqProj/results/firstPassAnalysis/K70MD_DEseq2.rda")
##===================================
##  K70 vs K70P10 (human)
##===================================
                              ##  exclude K70P10_3
#dds <- DESeqDataSetFromMatrix(countData=cleanMatrixData(d.cleaned [,c(1:4,8:10)]), 
# colData=metaData [c(1:3,7:9),] , 

dds <- DESeqDataSetFromMatrix(countData=cleanMatrixData(d.cleaned [,c(1:4,8:9)]), 
                              colData=metaData [c(1:3,7:8),] , 
                              design=~treatment, tidy = TRUE)
dds <- DESeq(dds)
res <- results(dds)
dim(res[(which(res$padj <= 0.05 )),])

##=========================================
##  examine the p-value distribution
##=========================================
hist(res$pvalue, col = "lavender", main = "WT vs Deletion", xlab = "p-values")
library(fdrtool)
FDR.res <- fdrtool(res$stat, statistic= "normal", plot = T)

FDR.res$param[1, "sd"]
res[,"padj"]  <- p.adjust(FDR.res$pval, method = "BH")


hist(FDR.res$pval, col = "royalblue4", 
     main = "WT vs Deletion, correct null model", xlab = "CORRECTED p-values")


table(res[,"padj"] < 0.1)
plotMA(res)
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

DEGs.P10 <- rownames(subset(res, padj<.05 & abs(log2FoldChange)>2))

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



