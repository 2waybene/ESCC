source("x:/R-project/customPackages/plotTools.R")
source("x:/R-project/customPackages/dataManipTools.R")
library(gplots)
library(readxl)
library(DESeq2)
library(ggplot2)
library(BiocStyle)
library(rmarkdown)
library(geneplotter)
library(plyr)
library(LSD)
library(RColorBrewer)
library(stringr)
library(topGO)
library(genefilter)
library(biomaRt)
library(dplyr)
library(EDASeq)
library(fdrtool)
library(org.Mm.eg.db)

#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("datatable")
#BiocManager::install("topGO")
#BiocManager::install("EDASeq")
#BiocManager::install("BioStyle")  ## do not work
#BiocInstaller::biocLite('BiocStyle') ## works


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