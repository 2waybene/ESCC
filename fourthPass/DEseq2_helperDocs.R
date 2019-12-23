
##  credit: https://support.bioconductor.org/p/79950/
library("BiocParallel")
register(MulticoreParam(2))
library("DESeq2")

countTable = read.table("Treatment_3_vs_Ctrl_3.txt",header=TRUE, row.names=1)
condition = factor(c("Ctrl_3","Treatment_3"))
libType = c("oneFactor","oneFactor")
condition <- relevel(condition, "Ctrl_3")
experiment_design=data.frame(row.names = colnames(countTable), condition, libType)



cds <- DESeqDataSetFromMatrix(countData = countTable, colData=experiment_design, design=~condition)


cds_DESeqED <- DESeq(cds,parallel = TRUE)
res <- results(cds_DESeqED,parallel = TRUE,alpha = 0.05, pAdjustMethod = "BH")
write.table(res,file = "Treatment_3_vs_Ctrl_3.differentialExpression.txt",row.names = TRUE,col.names = NA,append = FALSE, quote = FALSE, sep = "\t",eol = "\n", na = "NA", dec = ".")


normalizedReadCounts = counts(cds_DESeqED,normalized=TRUE)
write.table(normalizedReadCounts,file = "Treatment_3_vs_Ctrl_3.normalizedCounts.xls",row.names = TRUE,col.names = NA,append = FALSE, quote = FALSE, sep = "\t",eol = "\n", na = "NA", dec = ".")

