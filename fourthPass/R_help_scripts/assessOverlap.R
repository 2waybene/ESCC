##======================================================
#	File:   assessingImbalance
#	Author: Jianying Li
#	Initial coded for "liver" sample (FoleyWang.R) in May 2012
#	Date:   01/06/2014
##======================================================
source("x:/R-project/customPackages/dataManipTools.R")
source("x://R-project/customPackages/arraySeqTools.R")
source("x:/R-project/customPackages/plotTools.R")
source("X:/project2012/ASE/R-scripts/functions-ASE.R")
library(gplots)
library(Rlab)
library(affy)




##====================End of functions=========================================


samples <- c("Gadd45", "Hemox1", "IL6", "IL8")

file.dir <- "x:/project2015/DieselWes/LM_analysis"

corrMethod <- c("spearman", "pearson")
j = 1

LM_DEG_lists = c(list())

for (i in 1:length(samples))
{
	my.dir <- paste (file.dir, samples[i], corrMethod[j], sep = "/")	
	file.name <- paste ("LM_DEG_", samples[i], "_", corrMethod[j], ".txt", sep = "")
	file      <- paste (my.dir , file.name , sep = "/")
	dt <- read.table (file)
	LM_DEG_lists[[i]] <- dt$V1
}




VennDiagram <- draw.four.list.venndigram(LM_DEG_lists[[1]], LM_DEG_lists[[2]], LM_DEG_lists[[3]],  LM_DEG_lists[[4]], samples[[1]], samples[[2]], samples[[3]], samples[[4]])
grid.draw(VennDiagram$figure)



VennDiagram <- draw.three.list.venndigram(LM_DEG_lists[[1]], LM_DEG_lists[[2]], LM_DEG_lists[[3]],  samples[[1]], samples[[2]], samples[[3]])
grid.draw(VennDiagram$figure)
spearman.intersect.DEGs <- VennDiagram$union
pearson.intersect.DEGs <- VennDiagram$union



##====================Now, let's look at the overlap between selected sets================================


samples <- c("Gadd45", "Hemox1", "IL6", "IL8")

file.dir <- "x:/project2015/DieselWes/fc_analysis"
fn <- list.files (pattern = "_min3.txt", path = file.dir)
fn


LM_DEG_lists = c(list())

for (i in 1:4)
{
	my.dir <- file.dir	
	file.name <- paste ("LM_fc_up_", samples[i], "_min3.txt", sep = "")
	file      <- paste (my.dir , file.name , sep = "/")
	dt <- read.table (file)
	LM_DEG_lists[[i]] <- dt$V1
}



VennDiagram <- draw.four.list.venndigram(LM_DEG_lists[[1]], LM_DEG_lists[[2]], LM_DEG_lists[[3]],  LM_DEG_lists[[4]], samples[[1]], samples[[2]], samples[[3]], samples[[4]])
grid.draw(VennDiagram$figure)



LM_DEG_lists = c(list())

for (i in 1:4)
{
	my.dir <- file.dir	
	file.name <- paste ("LM_fc_down_", samples[i], "_min3.txt", sep = "")
	file      <- paste (my.dir , file.name , sep = "/")
	dt <- read.table (file)
	LM_DEG_lists[[i]] <- dt$V1
}



VennDiagram <- draw.four.list.venndigram(LM_DEG_lists[[1]], LM_DEG_lists[[2]], LM_DEG_lists[[3]],  LM_DEG_lists[[4]], samples[[1]], samples[[2]], samples[[3]], samples[[4]])
grid.draw(VennDiagram$figure)






