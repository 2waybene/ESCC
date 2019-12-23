##  credit: https://www.huber.embl.de/users/klaus/Teaching/DESeq2Predoc2014.html#inspection-and-correction-of-pvalues


##  Examine DESeq2Table
head(assay(DESeq2Table))
colSums(assay(DESeq2Table))
colData(DESeq2Table)
rowData(DESeq2Table)

mcols(rowData(DESeq2Table))

