##======================================================================
##  File  : DEseq2_analysis_v3.R
##  Author: Jianying Li
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
DESeq2Res = deseq2.results.K70MD$results

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


K70MD.DEGs2.adjp01.fc2 <- rownames(subset(res, padj<.01 & abs(log2FoldChange)>1))

length(K70MD.DEGs2.adjp01.fc2)
# 1727

##===================================
##  K70 vs K70P10 (human)
##===================================
##  exclude K70P10_3
#dds <- DESeqDataSetFromMatrix(countData=cleanMatrixData(d.cleaned [,c(1:4,8:10)]), 
# colData=metaData [c(1:3,7:9),] , 

#dds <- DESeqDataSetFromMatrix(countData=cleanMatrixData(d.cleaned [,c(1:4,8:9)]), 
#                              colData=metaData [c(1:3,7:8),] , 
#                              design=~treatment, tidy = TRUE)

#dds <- DESeq(dds)
#res <- results(dds)


#deseq2.results.K70P10 <- list (dat = cleanMatrixData(d.cleaned [,c(1:7)]), model = dds, results = res)
#save (deseq2.results.K70P10, file = "X:/project2019/RNAseqProj/results/thirdPassAnalysis/P10_DEseq2.rda")


load ("X:/project2019/RNAseqProj/results/thirdPassAnalysis/P10_DEseq2.rda")
res <- deseq2.results.K70P10$results
table(res$padj < 0.01)
#FALSE  TRUE 
#11069  9098 

##  third pass analysis

with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="Human RNAseq Volcano plot -- K70 vs K70P10", xlim=c(-3,3)))
# Add colored points: blue if padj<0.01, red if log2FC>1 and padj<0.05)
with(subset(res, padj<.05 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(res, padj<.05 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
mtext ("adjusted p-value < 0.05, log2FoldChange > 2")
dim(subset(res, padj<.05 & abs(log2FoldChange)>2))
#[1] 927    6

with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="Human RNAseq Volcano plot -- K70 vs K70P10", xlim=c(-3,3)))
# Add colored points: blue if padj<0.01, red if log2FC>1 and padj<0.05)
with(subset(res, padj<.01 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(res, padj<.01 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
mtext ("adjusted p-value < 0.01, log2FoldChange > 2")
dim(subset(res, padj<.01 & abs(log2FoldChange)>2))
#[1] 701    6


DESeq2Res = res

hist(DESeq2Res$pvalue, col = "lavender", main = "K70P10 vs K70", xlab = "p-values")


dim(DESeq2Res)
#[1] 25630     6
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
     main = "K70P10 vs K70 correct null model", xlab = "CORRECTED p-values")

DESeq2Res[,"padj"]  <- p.adjust(FDR.DESeq2Res$pval, method = "BH")


table(DESeq2Res[,"padj"] < 0.01)
#FALSE  TRUE 
#17290 2877


res = DESeq2Res
##  third pass analysis


with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="Human RNAseq Volcano plot -- K70 vs K70P10", xlim=c(-3,3)))
# Add colored points: blue if padj<0.01, red if log2FC>1 and padj<0.05)
with(subset(res, padj<.05 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(res, padj<.05 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
mtext ("adjusted p-value < 0.05, log2FoldChange > 2")
dim(subset(res, padj<.05 & abs(log2FoldChange)>2))
#[1] 263    6

with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="Human RNAseq Volcano plot -- K70 vs K70P10", xlim=c(-3,3)))
# Add colored points: blue if padj<0.01, red if log2FC>1 and padj<0.05)
with(subset(res, padj<.01 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(res, padj<.01 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
mtext ("adjusted p-value < 0.01, log2FoldChange > 2")
dim(subset(res, padj<.01 & abs(log2FoldChange)>2))
#[1] 221    6

with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="Human RNAseq Volcano plot -- K70 vs K70P10", xlim=c(-3,3)))
# Add colored points: blue if padj<0.01, red if log2FC>1 and padj<0.05)
with(subset(res, padj<.01 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(res, padj<.01 & abs(log2FoldChange)>1), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
mtext ("adjusted p-value < 0.01, log2FoldChange > 1")
dim(subset(res, padj<.01 & abs(log2FoldChange)>1))
#[1] 1990    6

K70P10.DEGs1.adjp01.fc2 <- rownames(subset(res, padj<.01 & abs(log2FoldChange)>1))

length(K70P10.DEGs1.adjp01.fc2)
# 1990

##================================
##  Working with the gene list
##================================

length(K70P10.DEGs1.adjp01.fc2)
# 1990

length(K70MD.DEGs2.adjp01.fc2)
# 1727

intersect(K70P10.DEGs1.adjp01.fc2, K70MD.DEGs2.adjp01.fc2)

library(readxl)

PAM50 <- read_excel("x:/project2019/RNAseqProj/doc/PAM50.xlsx")
pam50.gene.id <- PAM50$EntrezID

ESCC <- read_excel("x:/project2019/RNAseqProj/doc/Nrf2_genes_DHuang.xlsx")
ESCC_Nrf2 <- ESCC$`Gene Name`

mouse_ChipSeq_DEGs_ensembl <- read_excel("x:/project2019/RNAseqProj/doc/mosue_Nrf2_ChIP-seq_Nrf2K-vsKeap1KO_esophagus_JYL_w_human_homolog.xlsx")$"Human gene name"
length(mouse_ChipSeq_DEGs_ensembl)
mouse_ChipSeq_DEGs <- mouse_ChipSeq_DEGs_ensembl[!is.na(mouse_ChipSeq_DEGs_ensembl)]
length(mouse_ChipSeq_DEGs)

VennDiagram.3 <- draw.three.list.venndigram(K70P10.DEGs1.adjp01.fc2, K70MD.DEGs2.adjp01.fc2, ESCC_Nrf2, "P10", "MD", "ESCCNrf2")
grid.draw(VennDiagram.3$figure)


VennDiagram.3$unionList

VennDiagram.3 <- draw.three.list.venndigram(K70P10.DEGs1.adjp01.fc2, K70MD.DEGs2.adjp01.fc2, pam50.gene.id, "P10", "MD", "PAM50")
grid.draw(VennDiagram.3$figure)


DEGs.3.pass <- list (P10 = K70P10.DEGs1.adjp01.fc2, MD = K70MD.DEGs2.adjp01.fc2, PAM50 = pam50.gene.id , Nrf2 = ESCC_Nrf2)
save (DEGs.3.pass , file = "X:/project2019/RNAseqProj/results/thirdPassAnalysis/DEGs_3_lists.rda")

##===========================
##  Heatmaps
##===========================

load ("X:/project2019/RNAseqProj/results/firstPassAnalysis/K70ALL4_DEseq2.rda")
normalizedReadCountsAll = deseq2.results.K70.all4$norm.dat


data <- normalizedReadCountsAll[which(rownames(normalizedReadCountsAll) %in% intersect(K70P10.DEGs1.adjp01.fc2, K70MD.DEGs2.adjp01.fc2)), c(1:8)]
data <- normalizedReadCountsAll[which(rownames(normalizedReadCountsAll) %in% intersect(K70P10.DEGs1.adjp01.fc2, K70MD.DEGs2.adjp01.fc2)), ]
data <- normalizedReadCountsAll[which(rownames(normalizedReadCountsAll) %in%VennDiagram.3$unionList),c(1:8)]




heatmap.2(as.matrix(data),   scale="row", key=T, keysize=1.5,
          density.info="none", trace="none",cexCol=0.9, labRow=NA)


##===========================
##  Get DEGs lists
##===========================

load("X:/project2019/RNAseqProj/results/thirdPassAnalysis/DEGs_3_lists.rda")
write.table (DEGs.3.pass$P10, file = "X:/project2019/RNAseqProj/results/thirdPassAnalysis/P10_adjp01_fc2_DEGs.txt", sep = "\t" , row.names = TRUE, col.names = NA)
write.table (DEGs.3.pass$MD, file = "X:/project2019/RNAseqProj/results/thirdPassAnalysis/P10MD_adjp01_fc2_DEGs.txt", sep = "\t" , row.names = FALSE, col.names = FALSE)
write.table (DEGs.3.pass$MD, file = "X:/project2019/RNAseqProj/results/thirdPassAnalysis/K10MD_adjp01_fc2_DEGs.txt", sep = "\t" , row.names = FALSE)

K70P10.DEGs1.adjp01.fc2 = DEGs.3.pass$P10
K70MD.DEGs2.adjp01.fc2  = DEGs.3.pass$MD

common.P10.MD.DEGs.adjp01.fc2          = intersect(K70P10.DEGs1.adjp01.fc2,K70MD.DEGs2.adjp01.fc2  )
common.P10.ESCCNrf2.DEGs.adjp01.fc2    = intersect(K70P10.DEGs1.adjp01.fc2,  ESCC_Nrf2 )
common.K70MD.ESCCNrf2.DEGs.adjp01.fc2  = intersect(K70MD.DEGs2.adjp01.fc2,  ESCC_Nrf2 )
K70P10.DEGs1.adjp01.fc2.uniq   = K70P10.DEGs1.adjp01.fc2 [-which(K70P10.DEGs1.adjp01.fc2 %in% common.P10.MD.DEGs.adjp01.fc2)]
K70MD.DEGs2.adjp01.fc2.uniq   = K70MD.DEGs2.adjp01.fc2 [-which(K70MD.DEGs2.adjp01.fc2 %in% common.P10.MD.DEGs.adjp01.fc2)]

length(K70MD.DEGs2.adjp01.fc2.uniq)
#[1] 1196
length(K70P10.DEGs1.adjp01.fc2.uniq)
#[1] 1459
length(common.P10.MD.DEGs.adjp01.fc2 )
#531
length(K70P10.DEGs1.adjp01.fc2 )
#[1] 1990
length(K70MD.DEGs2.adjp01.fc2)
#[1] 1727
length(common.P10.ESCCNrf2.DEGs.adjp01.fc2)    
#43
length(common.K70MD.ESCCNrf2.DEGs.adjp01.fc2)  
#47

intersect(K70P10.DEGs1.adjp01.fc2, mouse_ChipSeq_DEGs)
# 209
intersect(K70MD.DEGs2.adjp01.fc2, mouse_ChipSeq_DEGs)
# 185

intersect(K70MD.DEGs2.adjp01.fc2, mouse_ChipSeq_DEGs)

VennDiagram.3 <- draw.three.list.venndigram(K70P10.DEGs1.adjp01.fc2, K70MD.DEGs2.adjp01.fc2, mouse_ChipSeq_DEGs, "P10", "MD", "MouseChIPSeq")
grid.draw(VennDiagram.3$figure)



names(K70MD.DEGs2.adjp01.fc2.uniq) <- "K70MD.DEGs2.adjp01.fc2.uniq"
write.table (as.vector(K70MD.DEGs2.adjp01.fc2.uniq), file = "X:/project2019/RNAseqProj/results/thirdPassAnalysis/unique_K10MD_adjp01_fc2_DEGs.txt", sep = "\t" , row.names = FALSE)
write.table (as.vector(K70P10.DEGs1.adjp01.fc2.uniq ), file = "X:/project2019/RNAseqProj/results/thirdPassAnalysis/unique_P10_adjp01_fc2_DEGs.txt", sep = "\t" , row.names = FALSE)
write.table (as.vector(common.P10.MD.DEGs.adjp01.fc2), file = "X:/project2019/RNAseqProj/results/thirdPassAnalysis/common_P10_K10MD_adjp01_fc2_DEGs.txt", sep = "\t" , row.names = FALSE)
write.table (as.vector(common.P10.ESCCNrf2.DEGs.adjp01.fc2), file = "X:/project2019/RNAseqProj/results/thirdPassAnalysis/common_P10_ESCCNrf2_adjp01_fc2_DEGs.txt", sep = "\t" , row.names = FALSE)
write.table (as.vector(common.K70MD.ESCCNrf2.DEGs.adjp01.fc2), file = "X:/project2019/RNAseqProj/results/thirdPassAnalysis/common_K10MD__ESCCNrf2_adjp01_fc2_DEGs.txt", sep = "\t" , row.names = FALSE)
write.table (as.vector(mouse_ChipSeq_DEGs), file = "X:/project2019/RNAseqProj/results/thirdPassAnalysis/mouse_ChIPSeq_human_homolog.txt", sep = "\t" , row.names = FALSE)
write.table (as.vector(VennDiagram.3$unionList), file = "X:/project2019/RNAseqProj/results/thirdPassAnalysis/P10_MD_mouse_ChIPSeq_human_homolog.txt", sep = "\t" , row.names = FALSE)

##==========================
##  PCA, 
##==========================




##==========================
##  Pathway analysis
##==========================



##==========================
##  Mouse orthologue
##==========================

## done


##==========================
##  AREs related
##==========================



##==========================
##  LINCS -- chemicals
##==========================











