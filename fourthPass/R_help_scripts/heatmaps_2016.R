##   http://mannheimiagoesprogramming.blogspot.com/2012/06/drawing-heatmaps-in-r-with-heatmap2.html

setwd("x:/learningR/learningHeatMap/datasets/")
data   = read.table ("testOut-TPM-TE.txt", header = TRUE, sep="\t")
View(data)
row.names(data) <- data[,1]
data <- data[,-1]


heatmap(as.matrix(data))  

heatmap(as.matrix(data), cexCol=0.7, labRow=NA)



library(gplots)
help(heatmap.2)

heatmap.2(as.matrix(data),   scale="row", key=T, keysize=1.5,
          density.info="none", trace="none",cexCol=0.9, labRow=NA)


heatmap.2(as.matrix(data), col=redgreen(75), scale="row", key=T, keysize=1.5,
          density.info="none", trace="none",cexCol=0.9, labRow=NA)


data = log2(data)
boxplot(data, cex.axis=0.8, las=2, main="Original distribution of data",
        ylab="Log2(Intensity)")  # Draw a boxplot.



# Normalization
library(preprocessCore)
data2 = normalize.quantiles(as.matrix(data))   # A quantile normalization. 


# Copy the row and column names from data to data2:
rownames(data2) = rownames(data)
colnames(data2) = colnames(data)


boxplot(data2, cex.axis=0.8, las=2, main="Distribution after normalization",
        ylab="Log2(Intensity)") 



library(limma)

design = cbind(Cell_A = c(1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0), # First 3 columns->Cell_A
               Cell_B = c(0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
               Cell_C = c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
               Cell_D = c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0),
               Cell_E = c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1)) # Last 3 columns->Cell_C
fit = lmFit(data2, design=design) # Fit the original matrix to the above design.
# We want to compare A vs. B, A vs. C and B vs. C
contrastsMatrix = makeContrasts("Cell_A-Cell_B","Cell_A-Cell_C","Cell_A-Cell_D","Cell_A-Cell_E",
                                "Cell_B-Cell_C","Cell_B-Cell_D","Cell_B-Cell_E",
                                "Cell_C-Cell_D","Cell_C-Cell_E","Cell_D-Cell_E",
                                levels = design) 


fit2 = contrasts.fit(fit, contrasts = contrastsMatrix) # Making the comparisons.
fit2 = eBayes(fit2) # Moderating the t-tetst by eBayes method.

str(fit2)

data3 = as.matrix(fit)
data3[is.na(data3)] <- 1



Label = c(rep("purple",250),rep("orange",250),rep("darkgreen",250),
          rep("brown",323))

data3 = data

heatmap.2(data3, col=redblue(256), dendrogram="both",
          scale="row", key=T, keysize=0.5, density.info="none",
          trace="none",cexCol=1.2, labRow=NA, RowSideColors=Label,
          lmat=rbind(c(5,0,4,0),c(3,1,2,0)), lhei=c(2.0,5.0),
          lwid=c(1.5,0.2,2.5,2.5))

x <- data3

heatmap.2((as.matrix(x)), col=redgreen(75), scale="row", key=T, keysize=1.5,
          density.info="none", trace="none",cexCol=0.9, labRow=NA)

heatmap.2((as.matrix(x)), col=heat.colors(100), scale="row", key=T, keysize=1.5,
          density.info="none", trace="none",cexCol=0.9, labRow=NA)


heatmap.2((as.matrix(x)), col=cm.colors(100), scale="row", key=T, keysize=1.5,
          density.info="none", trace="none",cexCol=0.9, labRow=NA)



heatmap.2((as.matrix(x)), col=terrain.colors(100), scale="row", key=T, keysize=1.5,
          density.info="none", trace="none",cexCol=0.9, labRow=NA)




##=============Okay upto here=============================
##
##  Now, need to see whether the following works
##  Need to create FoldCh and Pval two datamatrix
##
##========================================================
cmps <- c("Cell_A-Cell_B","Cell_A-Cell_C","Cell_A-Cell_D","Cell_A-Cell_E",
                     "Cell_B-Cell_C","Cell_B-Cell_D","Cell_B-Cell_E",
                     "Cell_C-Cell_D","Cell_C-Cell_E","Cell_D-Cell_E")
mirs <-  rownames(aa)



FoldCh <- matrix (nrow = dim(data3)[1], ncol=length(cmps))
Pval <- matrix (nrow = dim(data3)[1], ncol=length(cmps))



index = 1;
pcut = 1; 
fcut = 0

for (cmp in cmps)
{
#cmp <- cmps[1]


  res = decideTests(fit2,p.value=pcut,adjust.method="BH",lfc=log2(fcut))
  res = as.matrix(res)
  af.ids = names(res[abs(res[,cmp])==1,cmp])
  aa = topTable(fit2,coef=cmp,number=Inf,adjust="BH");
  FoldCh[,index] <- aa$logFC
  Pval [,index] <- aa$P.Value
  index = index + 1
}

rownames(FoldCh) =  mirs 
colnames(FoldCh) =  cmps
rownames(Pval) =  mirs 
colnames(Pval) =  cmps

pval.nas <- which(is.na(Pval[,1]))
Pval   <- Pval  [-pval.nas,]
FoldCh <- FoldCh[-pval.nas,]

##=====================================
##  simualted pvalue -- significant
##  -1: non-significant
##   0: significant
##   1: very-significant
##====================================
PVal   <- matrix (sample( rep(c(-1,0,1),100),100), nrow=10)
FoldCh <- FoldCh[c(1:10),]


#FoldCh = read.table("Fold_changes.txt",header=T,stringsAsFactors=F)
#Pval   = read.table("Pvalues.txt",header=T,stringsAsFactors=F)

heatmap.2(as.matrix(FoldCh),dendrogram="col", cellnote=as.matrix(PVal),
          notecol="black",col=redblue(256),scale="none",key=TRUE, keysize=1.5,
          density.info="none", trace="none", cexRow=0.7,cexCol=1.2)


##=====================================
##  pheatmap
##  pretty heatmap package
##====================================


library(pheatmap)   
library(gplots)
if (nrow(data) > 100) stop("Too many rows for heatmap, who can read?!")
fontsize_row = 10 - nrow(data) / 15
fontsize = 2
fontsize_row = 1
pheatmap(data, col=greenred(256), main="My Heatmap", cluster_cols=F, 
         fontsize_row=fontsize_row, border_color=NA)



