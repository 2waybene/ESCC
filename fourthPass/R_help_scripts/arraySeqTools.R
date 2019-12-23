##===============================================================================================================
# File name: arraySeqTools.R
# Author: Jianying Li
# History: initially coded by Joel Parker as arrayTools.R
#	     Modified by Jianying Li (Since December 2009 for specialized purposes
#	     Modified "readarray" method to accept just single array
##===============================================================================================================
### TO DO
#	annotate this file with output
# error when levels are 0,1
# write help for myHeatmap, qqPval, and updated readarray
# incorporate dlda (from supclust) as a predictor
# incorporate rparallel library
# dwd based functions should return
# build install file
###

#		batchAdjust(x,batchID,method="median",pcNum=NULL)
#				x - data matrix
#				batchID -	single row from a matrix of sample labels that indicate batch, may be numeric or text
#				method	-	method of transformation, one of median, mean, dwd, or pc
#				pcNum		-	the principal component to be removed if method="pc" is specified

#		clancCV(x,classes,geneSetSizes=c(seq(1,20,by=1)),priors="equal",folds=10)
#				x				- data data matrix
#				classes	-	a list labels for use in coloring the points
#				geneSetSizes	-	a vector of integers the will specify the number of genes in each class
#				priors	-	'equal' assumes equal class priors, 'class' calculates them based on proportion in the data
#				folds		-	the number of folds of CV

#		clancPredict(x,classes,y,ngenes,priors="equal",std=TRUE)
#				x				- the data matrix of training samples
#				classes	-	a list labels for use in coloring the points
#				y				-	the data matrix of test samples
#				ngenes	-	the number of genes selected when training the model
#				priors	-	'equal' assumes equal class priors, 'class' calculates them based on proportion in the data
#				std			- when true, the training and testing samples are standardized to mean=0 and var=1

#		co.var <- (x)
#				x				- a data vector
#

#		collapseIDs(x,allids,method="mean")
#				x				- a data matrix
#				allids	- a list of non-unique ids that correspond to the order of the data matrix
#				method	-	the method for collapsing data, can be one of mean, median, or sd

#		dwd(x,y)
#			x	- data matrix
#			y	-	class labels
				
#		dwdadjust(x,y,dwddir)
#			x	- data matrix
#			y	-	class labels
#			dwddir	-	the dwd direction vector output from dwd

#		dwdCV(x,y,cvFolds=10,nGenes=50)
#				x				- the data matrix of training samples
#				classes	-	a list labels for use in coloring the points
#				y				-	the data matrix of test samples	
#				cvFolds	-	the number of folds of CV
#				nGenes	-	the number of genes selected when training the model

#		dwdPredict(x,classes,y,std=F)
#				x				- the data matrix of training samples
#				classes	-	a list labels for use in coloring the points
#				y				-	the data matrix of test samples	
#				std			- when true, the training and testing samples are standardized to mean=0 and var=1

#		edge(x,classes=NA,times,permutations=500)
#				x				- a data matrix
#				classes	-	a list labels for use in coloring the points
#				times		-	a quantitative variable indicating time
#				permutations	-	the number of permutations of randomized class labels

#		edgeOneClass(x,classes=NA,times,permutations=500)
#				x				- a data matrix
#				classes	-	a list labels for use in coloring the points, required so that subsets of samples can be easily identified
#				times		-	a quantitative variable indicating time
#				permutations	-	the number of permutations of randomized class labels

#		fcAndT(x,classes)
#				x				- a data matrix
#				classes	-	a list of encoded group labels that define the experimental groups and the test.
#									the encoding scheme is identical to the requirements of SAM in Excel

#		kmPlot(x,event,stime,varName="",ymin=0,lineColors=NA,overall=F)
#				x				- vector of categorical group labels for each sample
#				event		-	vector of the event status for each sample
#				stime		- vector of the time to event or censor for each sample
#				varName	-	main title on the plot
#				ymin		-	minimum survival probability on the y-axis
#				lineColors	-	the line colors for the groups
#				overall	-	whent true, plot a reference line that encompasses all samples

#		knnCV(x,classes,geneSetSizes=c(seq(1,20,by=1)),nk=3)
#				x				- data data matrix
#				classes	-	a list labels for use in coloring the points
#				geneSetSizes	-	a vector of integers the will specify the number of genes in each class
#				nk			- k, the number of samples in the trained model	
 
#		knnPredict(x, classes, y, ngenes, nk)
#				x				- the data matrix of training samples, or pre-calculated centroids
#				classes	-	a list labels for use in coloring the points
#				y				-	the data matrix of test samples
#				ngenes	-	the number of genes selected when training the model
#				nk			- k, the number of samples in the trained model	
 
#		medianCtr(x)
#				x	-	a data matrix

#		multiGSA(x,classes,gmtFile,permutations=500,method="maxmean",minsize=15,logVals=T)
#				x				- data data matrix
#				classes	-	a list of encoded group labels that define the experimental groups and the test.
#									the encoding scheme is identical to the requirements of SAM in Excel
#				gmtFile	-	the file name containing a relationship of gene identifiers to category identifiers
#									gmt files are available from mSigDB (part of GSEA) and the GSA home page
#				permutations	-	the number of permutations of randomized class labels
#				method 	-	the method for calculating the univariate statistic, see GSA documentation for possible methods
#				minsize	-	the minimun number of genes in a set, sets with fewer genes are excluded
#				logVals	-	when true, it assumes the values are on log2 scale, otherwise this transformation is made

#		multiSam(x,classes,permutations=100,logVals=T)
#				x				- a data matrix
#				classes	-	the list of encoded group labels that define the experimental groups and the test.
#									the encoding scheme is taken from the SAM for Excel implementation:
#										two class paired: unique integers for each class
#										two class unpaired: integers designate the pairs and sign designates the class									
#										multiclass: 
#				permutations	-	the number of permutations of randomized class labels
#				logVals	-	when true, it assumes the values are on log2 scale, otherwise this transformation is made

#		overlapSets(x,y)
#				x	-	a data matrix
#				y	-	a data matrix

#		overlapVectors(x,y) #Added by Jianying Li for handling single vector case
#				x	-	a data frame (with name)
#				y	-	a data frame (with name)



#		pamCV(x,classes,folds=10,deltas=NULL)
#				x				- data data matrix
#				classes	-	a list labels for use in coloring the points
#				folds		-	the number of folds of CV
#				deltas	- a list of delta thresholds to use during CV

#		pca(x,classes,startPC=1,stopPC=4,jitt=5,size=1,legendloc="topright",returnLoadings=F)		
#				x					-	data matrix
#				classes 	- a list labels for use in coloring the points
#				startPC		- the first pc to be displayed
#				stopPC		- the last pc to be displayed
#				jitt			- the amount of random jitter applied to the y dimension of the points in the density plots
#				size			- size of the text and points
#				legendloc	-	the location of the legend on the upper left pc
#				returnLoadings	-	when true, a matrix of pc loadings is returned

#		pcaAndGSA(x,gmtFile,startPC=1,stopPC=4,fdrCut=0.05,permutations=200)
#				x					- data matrix
#				gmtFile		-	the file name containing a relationship of gene identifiers to category identifiers
#										gmt files are available from mSigDB (part of GSEA) and the GSA home page
#				startPC		-	the first pc to be tested
#				stopPC		- the last pc to be tested
#				fdrCut		- the FDR at which to stop printing results
#				permutations -	the number of permutations for the GSA algorithm

#		pcaEA(x,classes,size=1,showLegend=T,legendloc="topright",mainStr="",startPC=1,stopPC=2,showNames=T)
#				x					-	data matrix
#				classes 	- a list labels for use in coloring the points
#				size			- size of the text and points
#				showLegend	-	when true, show the legend it plot
#				legendloc	-	the location of the legend on the upper left pc
#				mainStr		-	the main title of the plot
#				startPC		- the first pc to be displayed
#				stopPC		- the last pc to be displayed
#				showNames	-	plot the sample names instead of points		

#		plotCV(x,k=1,mainTitle="")
#				x				- a list of cv results output from one of pamCV, clancCV, knnCV, or sspCV
#				k				- the index of the cv results to plot if x contains more than one CV result
#				mainTitle	-	the main title of the plot

#		plotdwd(x,y,dwddir)
#			x	- data matrix
#			y	-	class labels
#			dwddir	-	the dwd direction vector output from dwd

#		readarray(fileName,hr=1,impute=T,method="median")
#				fileName - 	name of the tab-delimited text file to read and parse. The file should have one column
#										of variable (i.e. gene) identifiers
#				hr			 -	the number of header rows, including the top row of sample labels
#				impute	 -	when true, imputation is performed via KNN imputation
#				method	 -	the method for collapsing duplicate gene IDs, can be one of mean, median, or sd

#		standardize(x)
#				x	-	a data matrix

#		sspCV(x,classes,geneSetSizes=c(seq(1,20,by=1)),folds=10)
#				x				- data data matrix
#				classes	-	a list labels for use in coloring the points
#				geneSetSizes	-	a vector of integers the will specify the number of genes in each class
#				priors	-	'equal' assumes equal class priors, 'class' calculates them based on proportion in the data
#				folds		-	the number of folds of CV

#		sspPredict(x,classes="",y,nGenes="",priors="equal",std=F,distm="euclidean",centroids=F)
#				x				- the data matrix of training samples, or pre-calculated centroids
#				classes	-	a list labels for use in coloring the points
#				y				-	the data matrix of test samples
#				nGenes	-	the number of genes selected when training the model
#				priors	-	'equal' assumes equal class priors, 'class' calculates them based on proportion in the data
#				std			- when true, the training and testing samples are standardized to mean=0 and var=1
#				distm 	-	the distance metric for determining the nearest centroid, can be one of euclidean, pearson, or spearman
#				centroids	-	when true, it is assumed that x consists of pre-calculated centroids

#		calcModules(x,gmtFile)
#				x				- data data matrix
#				gmtFile	-	the file name containing a relationship of gene identifiers to category identifiers
#									gmt files are available from mSigDB (part of GSEA) and the GSA home page

#		SWISS(data,groups,option,center) originally written by Chris Cabanski
#				data 		- (d x n) matrix with d genes and n samples
#				groups 	- (n x 1) vector with labels (1,2,3,...,k) for each group 
#							(i.e., cancer type)
#				option 	- option for centering
#						1 = usual sum of squares decomposition (uses overall mean)
#						2 = uses mean of the group means (this weighs all groups equally)
#				center 	- (d x 1) user input for centertering (for option = 3)
#						- if option = 1 or 2, then do not put anything
	
#		icc.2.vecs (seq.vector, ma.vector)
#				seq.vector    - a data vector to test
#		 		ma.vector     - a data vector to compare to
#				return: list of ICC, Rho, and data frame with two columns
 
#		subset.by.column (dm, columnHeaders)
#				dm		      - a data frame (nxm)
#		 		columnHeaders     - a string vector of column headers
#				return: a subset of dataframe having matched columnHeaders

###########################

#		dwdCV.inner(x,y,folds) INTERNAL
#		knnCV.internal(dm, tclass, cvIdx, active, nk)	INTERNAL
#		removepc(x,pc)	INTERNAL
#		sdfilter(x,sdt=1)	INTERNAL
#		chooseSAMtest(gr.labels) INTERNAL
#		cvSSP(dm, tclass, cvIdx, active) INTERNAL

###########################
library(samr)
library(GSA)
library(multtest)
#library(pamr)
#library(clanc) #train, build, test, no library
#source("C:/Users/li11/Documents/R-project/customPackages/clanc.R")
library(affy)
library(limma)
library(qvalue)
library(e1071)
library(class)
library(ctc)
library(supclust)
library(heatmap.plus)
library(irr)
#library(aroma.light)
#source("V:/current_projects/R_modules/dataManipTools.R")
#############################################################

qqPvals<-function(x,qs=c(0.95),xlab="Expected",ylab="Observed",main="QQ Plot",pch=19,col=1,cex=1){
	expected<- -log10(runif(sum(!is.na(x))))
	plot(sort(expected),sort(-log10(x)),xlab=xlab,ylab=ylab,main=main,pch=pch,cex=cex,col=col)
	abline(h=qs,lty=2)
}

cols <- function(lowi = "green", highi = "red", ncolors = 20) {
        low <- col2rgb(lowi)/255
        high <- col2rgb("black")/255
        col1 <- rgb(seq(low[1], high[1], len = ncolors), seq(low[2], 
            high[2], len = ncolors), seq(low[3], high[3], len = ncolors))
        low <- col2rgb("black")/255
        high <- col2rgb(highi)/255
        col2 <- rgb(seq(low[1], high[1], len = ncolors), seq(low[2], 
            high[2], len = ncolors), seq(low[3], high[3], len = ncolors))
        col<-c(col1[1:(ncolors-1)],col2)
        return(col)
}
  
myHeatmap<-function(x,t.colors=NA,fileName="cluster.cdt",linkage="complete",distance="pearson",contrast=2,returnSampleClust=F,rowNames=NA){
	
	temp<-hclust2treeview(x,method=distance,file=fileName,link=linkage,keep.hclust=T)
	gTree<-temp[[1]]
	sTree<-temp[[2]]
	
	imageVals<-x
	imageVals[x > contrast] <- contrast
	imageVals[x < -1 * contrast] <- -1 * contrast
	
	if(sum(is.na(t.colors))>0){
			heatmap(imageVals,Rowv=as.dendrogram(gTree),Colv=as.dendrogram(sTree),
							col=cols(),labCol=NA,scale="none",
							margins=c(1,7),labRow=rowNames)
	}else{
		if(length(t.colors)>dim(imageVals)[2]){
			heatmap.plus(imageVals,Rowv=as.dendrogram(gTree),Colv=as.dendrogram(sTree),
								col=cols(),labCol=NA,labRow=rowNames,scale="none",
								ColSideColors=t.colors, margins=c(1,7))
		}else{
			heatmap(imageVals,Rowv=as.dendrogram(gTree),Colv=as.dendrogram(sTree),
								col=cols(),labCol=NA,labRow=rowNames,scale="none",
								ColSideColors=as.vector(t(t.colors)), margins=c(1,7))
		}
	}
	if(returnSampleClust){
		return(sTree)
	}
}


batchAdjust<-function(x,batchID,method="median",pcNum=NULL){

	batchids <- levels(as.factor(as.vector(t(batchID))))
	allAnn<-dimnames(x)
	
	if(method=="median"){
		for(i in 1:length(batchids)){
			x[,batchids[i]==batchID] <- medianCtr(x[,batchids[i]==batchID])
		}	
	}
	if(method=="mean" | method=="anova"){
		for(i in 1:length(batchids)){
			x[,batchids[i]==batchID] <- t(scale(t(x[,batchids[i]==batchID]),scale=F))
		}	
	}
	if(method=="pca"){
		x<-removepc(x,pcNum)
	}

	dimnames(x)<-allAnn
	return(x)
}

multiSam<-function(x,classes,permutations=100,logVals=T){

	features<- dim(x)[1]
	sampleNames<- names(x)
	geneNames<-row.names(x)

	tempTable<-matrix()
	if(is.data.frame(classes)){
		rowsOfClasses <- dim(classes)[2]
		classes<-t(classes)
  	
		for(i in 1:rowsOfClasses){
  	
			gr.labels<-frmtClasses(classes[i,])
			dm <- x[,!(is.na(gr.labels))]
			samples<- dim(dm)[2]
		#stime<-classes[4,!(is.na(gr.labels))]
			names(dm) <- sampleNames[!(is.na(gr.labels))]
			gr.labels <- gr.labels[!(is.na(gr.labels))]
			whichtest<-chooseSAMTest(gr.labels)
			print(whichtest)
  	
			k<-nlevels(as.factor(gr.labels))
			zeros<-matrix(nrow=features,ncol=k)
			if(whichtest$t!="Quantitative"){
				if(whichtest$t=="Two class paired"){
					zeros<-matrix(nrow=features,ncol=2)
					zeros[,1]<-apply(dm[,gr.labels<0],1,sd)
					zeros[,2]<-apply(dm[,gr.labels>0],1,sd)
					minzero<-apply(zeros,1,min)
					dm[minzero<=0,]<-dm[minzero<=0,]+matrix(rnorm(samples,0,.0001),
								nrow=features,ncol=samples)[minzero<=0,]
				}else{
					for(j in 1:k){
						zeros[,j]<-apply(dm[,gr.labels==j],1,sd)
					}
					minzero<-apply(zeros,1,min)
					dm[minzero<=0,]<-dm[minzero<=0,]+matrix(rnorm(samples,0,.0001),
								nrow=features,ncol=samples)[minzero<=0,]
				}
			}
  	
			x.sam<-list(x=dm,y=gr.labels,geneid=geneNames,genenames=geneNames, logged2=logVals)
			x.sam.out<-samr(x.sam, resp.type=whichtest$t, nperms=permutations)
			
			#testtype<-"Survival"
			#x.sam<-list(x=dm,y=as.vector(t(stime)),geneid=geneNames,genenames=geneNames, logged2=logVals, censoring.status=gr.labels)
			#x.sam.out<-samr(x.sam, resp.type=testtype, nperms=permutations)
  	
			delta.table<-samr.compute.delta.table(x.sam.out)
  	
			siggenes.table<-samr.compute.siggenes.table(x.sam.out, 0, x.sam, delta.table,all.genes=T)
			temp <- rbind(siggenes.table$genes.up,siggenes.table$genes.lo)
			temp <- temp[sort.list(temp[,"Row"]),]
			if(i==1){
				tempTable <- cbind(temp[,"Gene ID"],temp[,"Gene Name"],temp[,whichtest$s],temp[,"q-value(%)"])
				dimnames(tempTable)[[2]] <- c("Gene ID","Gene Name",
						paste(rownames(classes)[i],whichtest$s),
						paste(rownames(classes)[i]," q-value(%)"))
			}else{
				tempTable <- cbind(tempTable,temp[,whichtest$s],temp[,"q-value(%)"])
				tempDim <- dim(tempTable)[2]
				dimnames(tempTable)[[2]][(tempDim-1):tempDim] <-
					c(paste(rownames(classes)[i],whichtest$score),paste(row.names(classes)[i]," q-value(%)"))
			}
		}	
	}else{
		rowsOfClasses <- 1
		gr.labels<-frmtClasses(classes)
		dm <- x[,!(is.na(gr.labels))]
		samples<- dim(dm)[2]
		#stime<-classes[4,!(is.na(gr.labels))]
		names(dm) <- sampleNames[!(is.na(gr.labels))]
		gr.labels <- gr.labels[!(is.na(gr.labels))]
		whichtest<-chooseSAMTest(gr.labels)
		print(whichtest)
  	
		k<-nlevels(as.factor(gr.labels))
		zeros<-matrix(nrow=features,ncol=k)
		if(whichtest$t!="Quantitative"){
			if(whichtest$t=="Two class paired"){
				zeros<-matrix(nrow=features,ncol=2)
				zeros[,1]<-apply(dm[,gr.labels<0],1,sd)
				zeros[,2]<-apply(dm[,gr.labels>0],1,sd)
				minzero<-apply(zeros,1,min)
				dm[minzero<=0,]<-dm[minzero<=0,]+matrix(rnorm(samples,0,.0001),
							nrow=features,ncol=samples)[minzero<=0,]
			}else{
				for(j in 1:k){
					zeros[,j]<-apply(dm[,gr.labels==j],1,sd)
				}
				minzero<-apply(zeros,1,min)
				dm[minzero<=0,]<-dm[minzero<=0,]+matrix(rnorm(samples,0,.0001),
							nrow=features,ncol=samples)[minzero<=0,]
			}
		}
  	
		x.sam<-list(x=dm,y=gr.labels,geneid=geneNames,genenames=geneNames, logged2=logVals)
		x.sam.out<-samr(x.sam, resp.type=whichtest$t, nperms=permutations)
		
		#testtype<-"Survival"
		#x.sam<-list(x=dm,y=as.vector(t(stime)),geneid=geneNames,genenames=geneNames, logged2=logVals, censoring.status=gr.labels)
		#x.sam.out<-samr(x.sam, resp.type=testtype, nperms=permutations)
  	
		delta.table<-samr.compute.delta.table(x.sam.out)
  	
		siggenes.table<-samr.compute.siggenes.table(x.sam.out, 0, x.sam, delta.table,all.genes=T)
		temp <- rbind(siggenes.table$genes.up,siggenes.table$genes.lo)
		temp <- temp[sort.list(temp[,"Row"]),]

		tempTable <- cbind(temp[,"Gene ID"],temp[,"Gene Name"],temp[,whichtest$s],temp[,"q-value(%)"])
		dimnames(tempTable)[[2]] <- c("Gene ID","Gene Name",
																		paste("SAM",whichtest$s),
																		paste("SAM"," q-value(%)"))
	}	

	return(tempTable)
}

standardize<-function(x){
	annAll<-dimnames(x)
	x<-scale(x)
	dimnames(x)<-annAll
	return(x)
}

pca<-function(x,classes,startPC=1,stopPC=4,jitt=5,size=1,legendloc="topright",returnLoadings=F){

	features<- dim(x)[1]
	samples<- dim(x)[2]
	sampleNames<- names(x)
	featureNames<-row.names(x)
	x<-apply(x,2,as.numeric)

	#principal components plots
	data.pca<-prcomp(as.matrix(x))

	# Proportion of total variance distributed over 10 first components:
	tmp<-data.pca$sdev[1:10]^2/sum(data.pca$sdev^2)

	gr.labels<-as.vector(classes)
	gr.labels.fac<-factor(classes,exclude="")
	nlabels<-nlevels(gr.labels.fac)
	legendLabels<-vector()
	for(k in 1:nlabels){
		group<-levels(gr.labels.fac)[k]
		legendLabels[k]<-group
		gr.labels[gr.labels.fac==group]<-k
	}
	gr.labels<-as.numeric(gr.labels)

# edit for more than 8 groups - change pch

	#plot 2pcs by each other
	pcs <- stopPC - startPC + 1
	par(mfrow=c(pcs,pcs))
	for( i in startPC:stopPC ){
		for(j in startPC:stopPC){

			#graphing parameters
			par(lab=c(3,4,3))
			par(mgp=c(.3,.5,.0))
			par(mai=c(.2,.2,.2,.2))
			par(tcl=c(-.2))
			par(xaxt="n",yaxt="n")

			if(j == i){
				str<-paste("PC",i,round(tmp[i],2),sep=" ")
				den<-density(data.pca$rotation[,i])
				y<-jitter(rep(max(den$y)/2,samples),jitt)
				plot(den,main=str,type="l",xlab="",ylab="Frequency")
				points(x=data.pca$rotation[,i],y=y,col=gr.labels,cex=size,pch=19)
				if(j==startPC){
					legend(legendloc,legend=legendLabels,col=seq(1,nlabels),pch=19,
									x.intersp=.3,yjust=.5,bty="n",cex=.8)
				}
			}else{
				strM<-paste("PC",i,"vs",j,sep=" ")
				strX<-paste("PC",i,sep=" ")
				strY<-paste("PC",j,sep=" ")
				plot(data.pca$rotation[,i],data.pca$rotation[,j], xlab=strX,pch=19,
					ylab=strY, col=gr.labels,cex=size)

			}
		}
	}
	if(returnLoadings){
		return(data.pca$rotation)
	}
}

medianCtr<-function(x){
	annAll <- dimnames(x)
	medians <- apply(x,1,median,na.rm=T)
	x <- t(scale(t(x),center=medians,scale=F))
	dimnames(x) <- annAll
	return(x)
}

removepc<-function(x,pc){
	s<-svd(x)
	D<-diag(s$d)
	for(i in 1:length(pc)){
		D[pc[i],pc[i]]<-0
	}
	dataTR<-s$u%*%D%*%t(s$v)
	return(dataTR)
}

readarray<-function(dataFile,designFile=NA,hr=1,impute=T,method="mean"){

	headerRows <- hr

	x<-read.table(dataFile,sep="\t",header=F,fill=T,stringsAsFactors=FALSE)

	if(headerRows==1){
			sampleNames<-as.vector(t(x[1,-1]))
			x<-x[-1,]
			classes<-NULL
			ids<-x[,1]
			xd<-x[,-1]
			xd<-as.matrix(xd) #modified by JYL for taking single array
			xd<-apply(xd,2,as.numeric)
			xd<-collapseIDs(xd,ids,method)	
	}else{
			sampleNames<-as.vector(t(x[1,-1]))
			x<-x[-1,]
			
			classes<-x[1:(headerRows-1),]
			dimnames(classes)[[1]]<-classes[,1]
			classes<-classes[,-1]
			classes[classes==""]<-NA
			classes<-t(classes)
			rownames(classes)<-sampleNames
			classes<-as.data.frame(classes)
						
			xd<-x[(-1:-(headerRows-1)),]
			ids<-as.vector(t(xd[,1]))
			xd<-xd[,-1]
			xd<-as.matrix(xd) #modified by JYL for taking single array
			xd<-apply(xd,2,as.numeric)
			xd<-collapseIDs(xd,ids,method)
	}
	
	features<- dim(xd)[1]
	samples<- dim(xd)[2]
	geneNames<-rownames(xd)
	xd<-apply(xd,2,as.numeric)
	rownames(xd)<-geneNames
	colnames(xd)<-sampleNames

	if(!is.na(designFile)){
		x<-read.table(designFile,sep="\t",header=T,row.names=1,fill=T,,stringsAsFactors=FALSE)
		xd<-xd[,sort.list(colnames(xd))]
		xd<-xd[,colnames(xd) %in% rownames(x)]
		x<-x[rownames(x) %in% colnames(xd),]
		x<-x[sort.list(rownames(x)),]
		classes<-as.data.frame(x)
	}
	
	if(sum(apply(xd,2,is.na))>0 & impute){
		library(impute)
		allAnn<-dimnames(xd)
		data.imputed<-impute.knn(as.matrix(xd))
		xd<-data.imputed[1:features,]
		dimnames(xd)<-allAnn
	}
	
	return(list(xd=xd, classes=classes, nfeatures=features, nsamples=samples, fnames=geneNames, snames=sampleNames))
}

overlapSets<-function(x,y){

	# subset the two lists to have a commonly ordered gene list
	x<-x[dimnames(x)[[1]] %in% dimnames(y)[[1]],]
	y<-y[dimnames(y)[[1]] %in% dimnames(x)[[1]],]

	#and sort such that thing are in the correct order
	x<-x[sort.list(row.names(x)),]
	y<-y[sort.list(row.names(y)),]

	return(list(x=x,y=y))
}


multiGSA<-function(x,classes,gmtFile,permutations=500,method="maxmean",minsize=15,logVals=T,whichtest=NULL){

	features<- dim(x)[1]
	samples<- dim(x)[2]
	sampleNames<- names(x)
	geneNames<-row.names(x)
	x<-apply(x,2,as.numeric)

	geneset.obj<- GSA.read.gmt(gmtFile)

	tempTable<-matrix()
	if(is.data.frame(classes)){
		rowsOfClasses <- dim(classes)[2]
		classes<-t(classes)

		for(i in 1:rowsOfClasses){
  	
			gr.labels<-frmtClasses(classes[i,])
			dm <- x[,!(is.na(gr.labels))]
			names(dm) <- sampleNames[!(is.na(gr.labels))]
			gr.labels<-gr.labels[!(is.na(gr.labels))]
			if(length(whichtest$t)==0){
				whichtest<-chooseSAMTest(gr.labels)
			}
			
			k<-nlevels(as.factor(gr.labels))
			zeros<-matrix(nrow=features,ncol=k)
			if(whichtest$t=="Two class paired"){
				zeros<-matrix(nrow=features,ncol=2)
				zeros[,1]<-apply(dm[,gr.labels<0],1,sd)
				zeros[,2]<-apply(dm[,gr.labels>0],1,sd)
				minzero<-apply(zeros,1,min)
				dm[minzero<=0,]<-dm[minzero<=0,]+matrix(rnorm(samples,0,.0001),
							nrow=features,ncol=samples)[minzero<=0,]
			}else{
				for(j in 1:k){
					zeros[,j]<-apply(dm[,gr.labels==j],1,sd)
				}
				minzero<-apply(zeros,1,min)
				dm[minzero<=0,]<-matrix(rnorm(samples,0,.0001),nrow=features,ncol=dim(dm)[2])[minzero<=0,]
			}
  	
			GSA.obj<-GSA(dm,gr.labels, genenames=geneNames, method=method, minsize=minsize,s0.perc=-1,
										genesets=geneset.obj$genesets,  resp.type=whichtest$t, nperms=permutations)
			temp<-rbind(GSA.listsets(GSA.obj, geneset.names=geneset.obj$geneset.names,FDRcut=100,maxchar=255)$positive,
						GSA.listsets(GSA.obj, geneset.names=geneset.obj$geneset.names,FDRcut=100,maxchar=255)$negative)
  	
			temp <- temp[sort.list(temp[,"Gene_set"]),]
			if(i==1){
			tempTable <- cbind(temp[,"Gene_set"],temp[,"Gene_set_name"],temp[,"Score"],temp[,"p-value"],temp[,"FDR"])
			dimnames(tempTable)[[2]] <- c("Set ID","Set Name",
					paste(row.names(classes)[i],"Score"),
					paste(row.names(classes)[i],"p-value"),
					paste(row.names(classes)[i],"FDR"))
			}else{
				tempTable <- cbind(tempTable,temp[,"Score"],temp[,"p-value"],temp[,"FDR"])
				tempDim <- dim(tempTable)[2]
				dimnames(tempTable)[[2]][(tempDim-2):tempDim] <-
					c(paste(row.names(classes)[i],"Score"),
						paste(row.names(classes)[i],"p-value"),
						paste(row.names(classes)[i],"FDR"))
			}
		}
	}else{
			gr.labels<-frmtClasses(classes)
			dm <- x[,!(is.na(gr.labels))]
			names(dm) <- sampleNames[!(is.na(gr.labels))]
			gr.labels<-gr.labels[!(is.na(gr.labels))]
			if(length(whichtest$t)==0){
				whichtest<-chooseSAMTest(gr.labels)
			}
			
			k<-nlevels(as.factor(gr.labels))
			zeros<-matrix(nrow=features,ncol=k)
			if(whichtest$t=="Two class paired"){
				zeros<-matrix(nrow=features,ncol=2)
				zeros[,1]<-apply(dm[,gr.labels<0],1,sd)
				zeros[,2]<-apply(dm[,gr.labels>0],1,sd)
				minzero<-apply(zeros,1,min)
				dm[minzero<=0,]<-dm[minzero<=0,]+matrix(rnorm(samples,0,.0001),
							nrow=features,ncol=samples)[minzero<=0,]
			}else{
				for(j in 1:k){
					zeros[,j]<-apply(dm[,gr.labels==j],1,sd)
				}
				minzero<-apply(zeros,1,min)
				dm[minzero<=0,]<-matrix(rnorm(samples,0,.0001),nrow=features,ncol=dim(dm)[2])[minzero<=0,]
			}
  	
			GSA.obj<-GSA(dm,gr.labels, genenames=geneNames, method=method, minsize=minsize,s0.perc=-1,
										genesets=geneset.obj$genesets,  resp.type=whichtest$t, nperms=permutations)
			temp<-rbind(GSA.listsets(GSA.obj, geneset.names=geneset.obj$geneset.names,FDRcut=100,maxchar=255)$positive,
						GSA.listsets(GSA.obj, geneset.names=geneset.obj$geneset.names,FDRcut=100,maxchar=255)$negative)
  	
			temp <- temp[sort.list(temp[,"Gene_set"]),]
			tempTable <- cbind(temp[,"Gene_set"],temp[,"Gene_set_name"],temp[,"Score"],temp[,"p-value"],temp[,"FDR"])
			dimnames(tempTable)[[2]] <- c("Set ID","Set Name",
																		paste("GSA Score"),
																		paste("GSA p-value"),
																		paste("GSA FDR"))
	}

	return(tempTable)
}


sdfilter<-function(x,sdt=1){
	sds <- apply(x,1,sd,na.rm=T)
	return (x[sds>sdt,])
}

chooseSAMTest<-function(gr.labels){
	# tests not implemented: 	"Quantitative","Two class paired timecourse",
	# 				"Survival", "Pattern discovery"
	if(length(levels(as.factor(gr.labels))) == 1){
			test<-"One class"
			score<-"Numerator(r)"
	}
	if(length(levels(as.factor(gr.labels))) == 2){
			test<-"Two class unpaired"
			score<-"Fold Change"
	}
	if(length(levels(as.factor(gr.labels))) > 2){
			test<-"Multiclass"
			score<-"Numerator(r)"
	}
	if(length(levels(as.factor(gr.labels))) == length(gr.labels)){
			test<-"Quantitative"
			score<-"Numerator(r)"
	}
	if(sum(gr.labels<0,na.rm=T) > 0){
			test<-"Two class paired"
			score<-"Fold Change"
	}
	if(sum(gr.labels=="1Time1Start",na.rm=T) > 0){
			test<-"One class timecourse"
			score<-"Numerator(r)"
	}
	if(sum(gr.labels=="2Time1Start",na.rm=T) > 0){
			test<-"Two class unpaired timecourse"
			score<-"Numerator(r)"
	}
	return(list(t=test,s=score))
}

frmtClasses<-function(classes){
		if(sum(classes=="1Time1Start",na.rm=T) == 0){
			gr.labels <- as.vector(as.numeric(t(classes)))
		}else{
			gr.labels <- as.vector(t(classes))
			gr.labels[gr.labels==""] <- NA
		}
		return(gr.labels)
}


fcAndT<-function(x,classes){

  features<- dim(x)[1]
  samples<- dim(x)[2]
  sampleNames<- names(x)
  geneNames<-row.names(x)
  x<-apply(x,2,as.numeric)

	tempTable<-matrix()
	if(is.data.frame(classes)){
		rowsOfClasses <- dim(classes)[2]
		classes<-t(classes)

    for(i in 1:rowsOfClasses){
     gr.labels <- as.vector(as.numeric(t(classes[i,])))
     dm <- x[,!(is.na(gr.labels))]
     names(dm) <- sampleNames[!(is.na(gr.labels))]
     gr.labels <- gr.labels[!(is.na(gr.labels))]

     if(length(levels(as.factor(gr.labels))) == 2){
         test<-"t"
         score<-"Log Ratio"
         if(min(gr.labels)==1){
              gr.labels[gr.labels==1]<-0
              gr.labels[gr.labels==2]<-1
         }
     }
     if(length(levels(as.factor(gr.labels))) > 2){
         test<-"f"
         score<-"Numerator(r)"
         if(min(gr.labels)==1){
             tmp<-max(gr.labels)
             for(j in 1:tmp){
                     gr.labels[gr.labels==j]<-j-1
             }
         }
     }
     if(sum(gr.labels<0,na.rm=T) > 0){
         test<-"pairt"
         score<-"Log Ratio"
         tmp<-apply(cbind(gr.labels,(-1*gr.labels)),1,max)
         gr.labels<-gr.labels[sort.list(tmp)]
         gr.labels[gr.labels<0]<-0
         gr.labels[gr.labels>0]<-1
         dm <- dm[,sort.list(tmp)]
         names(dm) <- names(dm)[sort.list(tmp)]
     }
     
  	k<-nlevels(as.factor(gr.labels))
		zeros<-matrix(nrow=features,ncol=k)

				for(j in 1:k){
					zeros[,j]<-apply(dm[,gr.labels==(j-1)],1,sd)
				}
				minzero<-apply(zeros,1,min)
				dm[minzero<=0,]<-dm[minzero<=0,]+matrix(rnorm(dim(dm)[2],0,.0001),
							nrow=features,ncol=dim(dm)[2])[minzero<=0,]
		
     result<-mt.teststat(dm,gr.labels,test=test)
     if(test=="f"){
	     result<-1-pf(result,tmp-1,dim(dm)[2]-tmp)
     }else{
     	result<-2*(1-pt(abs(result),dim(dm)[2]-1))
     }
     result.score<-mt.teststat.num.denum(dm,gr.labels,test=test)
     fdr<-qvalue(result)$qvalues
     if(i==1){
             tempTable <- cbind(geneNames,result.score$teststat.num,result,fdr)
             dimnames(tempTable)[[2]] <- c("Gene Name",
                             paste(row.names(classes)[i],score),
                             paste(row.names(classes)[i],"p-value"),
                             paste(row.names(classes)[i],"q-value"))
     }else{
             tempTable <- cbind(tempTable,result.score$teststat.num,result,fdr)
             tempDim <- dim(tempTable)[2]
             dimnames(tempTable)[[2]][(tempDim-2):tempDim] <-
                     c(paste(row.names(classes)[i],score),
                     		paste(row.names(classes)[i],"p-value"),
                     		paste(row.names(classes)[i],"q-value"))
     }
   }
  }else{
  	 gr.labels <- as.vector(classes)
     dm <- x[,!(is.na(gr.labels))]
     names(dm) <- sampleNames[!(is.na(gr.labels))]
     gr.labels <- gr.labels[!(is.na(gr.labels))]

     if(length(levels(as.factor(gr.labels))) == 2){
         test<-"t"
         score<-"Log Ratio"
         if(min(gr.labels)==1){
              gr.labels[gr.labels==1]<-0
              gr.labels[gr.labels==2]<-1
         }
     }
     if(length(levels(as.factor(gr.labels))) > 2){
         test<-"f"
         score<-"Numerator(r)"
         if(min(gr.labels)==1){
             tmp<-max(gr.labels)
             for(j in 1:tmp){
                     gr.labels[gr.labels==j]<-j-1
             }
         }
     }
     if(sum(gr.labels<0,na.rm=T) > 0){
         test<-"pairt"
         score<-"Log Ratio"
         tmp<-apply(cbind(gr.labels,(-1*gr.labels)),1,max)
         gr.labels<-gr.labels[sort.list(tmp)]
         gr.labels[gr.labels<0]<-0
         gr.labels[gr.labels>0]<-1
         dm <- dm[,sort.list(tmp)]
         names(dm) <- names(dm)[sort.list(tmp)]
     }
     
  	k<-nlevels(as.factor(gr.labels))
		zeros<-matrix(nrow=features,ncol=k)

				for(j in 1:k){
					zeros[,j]<-apply(dm[,gr.labels==(j-1)],1,sd)
				}
				minzero<-apply(zeros,1,min)
				dm[minzero<=0,]<-dm[minzero<=0,]+matrix(rnorm(dim(dm)[2],0,.0001),
							nrow=features,ncol=dim(dm)[2])[minzero<=0,]
		
     result<-mt.teststat(dm,gr.labels,test=test)
     if(test=="f"){
	     result<-1-pf(result,tmp-1,dim(dm)[2]-tmp)
     }else{
     	result<-2*(1-pt(abs(result),dim(dm)[2]-1))
     }
     result.score<-mt.teststat.num.denum(dm,gr.labels,test=test)
     fdr<-qvalue(result)$qvalues
     tempTable <- cbind(geneNames,result.score$teststat.num,result,fdr)
		dimnames(tempTable)[[2]] <- c("Gene ID","Gene Name",
																		paste("t-test",whichtest$s),
																		paste("t-test"," q-value(%)"))
  }
   return(tempTable)
}

calcModules<-function(x, gmtFile, method="mean", scaleGenes=F){

	geneset.obj<- GSA.read.gmt(gmtFile)
	genenames<-row.names(x)
	np=length(geneset.obj$genesets)

	val=matrix(NA,nrow=np,ncol=ncol(x))
	if(scaleGenes){xs=t(scale(t(x),center=T,scale=T))}else{xs<-x}
	for(i in 1:np){
	 	gene.set=which(genenames %in% geneset.obj$genesets[[i]])
	 	gene.set=gene.set[!is.na(gene.set)]
		if(length(gene.set)>1){
		 	if(method=="mean"){
			 	val[i,]=colSums(xs[gene.set,,drop=F])/length(gene.set)
			}
			if(method=="median"){
				xs<-x
				val[i,]=t(apply(xs[gene.set,,drop=F],2,median))
			}
			if(method=="pca"){
				xs<-x
				y<-prcomp(as.matrix(xs[gene.set,,drop=F]))
				val[i,]=y$rotation[,1]
			}
		}
	}
	dimnames(val)<-list(geneset.obj$geneset.names,dimnames(x)[[2]])
	return(val)
}

pcaAndGSA<-function(x,gmtFile,startPC=1,stopPC=4,fdrCut=0.05,permutations=200){

# Does "loadings" retain the same column order as x?

	loadings<-pca(x,matrix(1,nrow=1,ncol=dim(x)[2]),startPC,stopPC,returnLoadings=T)
	sigSet<-list()
	for(i in startPC:stopPC){
		gsaOut<-multiGSA(x,matrix(loadings[,i],nrow=1,ncol=dim(x)[2]),gmtFile,permutations=permutations,whichtest=list(t="Quantitative",s="Correlation"))
		sigSet[[i]]<-gsaOut[gsaOut[,5]<=fdrCut,2]
	}
	return(sigSet)
}
		
pamCV<-function(x,classes,folds=10,deltas=NULL){

	dataMatrix<-x
	features<- dim(dataMatrix)[1]
	samples<- dim(dataMatrix)[2]
	sampleNames<- names(dataMatrix)
	geneNames<-row.names(dataMatrix)
	
	outTables<-list()
	
	tclass <- as.vector(classes)
	dm <- dataMatrix[,!(is.na(tclass))]
	tclass <- tclass[!(is.na(tclass))]
	
	nclasses<-nlevels(as.factor(tclass))
	classLevels<-levels(as.factor(tclass))
	for(j in 1:nclasses){
		tclass[tclass==classLevels[j]] <- j
	}
	tclass<-as.numeric(tclass)
	
	trainData <- list(x=dm,y=tclass, geneid=geneNames, genenames=geneNames)
	mytrain <- pamr.train(trainData,threshold=deltas)
	#new.scales <- pamr.adaptthresh(mytrain)
	#mytrain2 <- pamr.train(trainData, threshold.scale=new.scales)
	cvResults <- pamr.cv(mytrain, trainData, folds=balancedFolds(tclass, folds))
	
	minE<-1
	for(j in 1:length(cvResults$threshold)){
		if(cvResults$error[j]<minE){
			minE<-cvResults$error[j]
			thr <- cvResults$threshold[j]
		}
	}
	
	cvClassErrors<-matrix(nrow=length(cvResults$error),ncol=nclasses)
	for(j in 1:length(cvResults$error)){
		for(k in 1:nclasses){
			cvClassErrors[j,k]<-1-(sum(cvResults$yhat[tclass==k,j]==tclass[tclass==k])/length(tclass[tclass==k]))
		}
	}
	outTable<-cbind(cvResults$threshold,cvResults$size,cvResults$error,cvClassErrors)
	dimnames(outTable)<-list(NULL,c("Delta","nGenes","Overall Error",paste(classLevels,"Error")))
	
	return(outTable)
}


clancCV<-function(x,classes,geneSetSizes=c(seq(1,20,by=1)),priors="equal",folds=10){

#Here, geneSetSizes can't be bigger so that geneSetSizes*5 > dim(dataMatrix)[1]
#Fixed by Jianying Li, 10/22/2012


	dataMatrix<-x
	features<- dim(dataMatrix)[1]
	samples<- dim(dataMatrix)[2]
	sampleNames<- names(dataMatrix)
	geneNames<-row.names(dataMatrix)

	if (geneSetSizes*5 > features || length(geneSetSizes)*5 > features)
	{
         geneSetSizes = floor(features/5)
	}

	# format the dependent variable
	tclass <- as.vector(classes)
	nclasses<-nlevels(as.factor(tclass))
	classLevels<-levels(as.factor(tclass))
	for(j in 1:nclasses){
		tclass[tclass==classLevels[j]] <- j
	}
	tclass<-as.numeric(tclass)
	dm <- dataMatrix[,!(is.na(tclass))]
	tclass <- tclass[!(is.na(tclass))]
		
	# run CV

##===============This is NOT working with local clanc.R===========================================================#
#	cvResults<-cvClanc(dm, tclass, cvIdx = balancedFolds(tclass, folds), prior = priors, active = geneSetSizes) #
#	Use the following function call instead												#
##================================================================================================================#

	cvResults<-cvClanc(dm, tclass, prior = priors, active = geneSetSizes)
	
	outTable<-cbind(geneSetSizes,geneSetSizes*nclasses,cvResults$overallErrors,cvResults$classErrors)
	dimnames(outTable)<-list(NULL,c("nGenes per class","nGenes","Overall Error",paste(classLevels,"Error")))
		
	return(outTable)
}

plotCV<-function(x,mainTitle=""){
	nclasses<-dim(x)[2]-3
	print(nclasses)
	plot(x[,2],x[,3],col=1,lty=1,lwd=3,ylim=c(0,.5),type="l",xlab="Features",ylab="CV Error",main=mainTitle)
	for(j in 1:nclasses){
		points(x[,2],x[,(j+3)],col=j,type="l",lty=1)
	}
	legend("topright",legend=dimnames(x)[[2]][3:(nclasses+3)],col=c(1,seq(1,nclasses)),lty=rep(1,(nclasses+1)),lwd=c(3,rep(1,nclasses)),bty="n")
}

clancPredict<-function(x,classes,y,ngenes,priors="equal",std=TRUE){
	activeGenes <- ngenes

	dataMatrix<-x
	sampleNames<- dimnames(x)[[2]]
	featureNames<- dimnames(x)[[1]]
	
	#parse the test file - same as train file but no rows of classes
	tdataMatrix<-y
	tsampleNames<- dimnames(y)[[2]]
	tfeatureNames<- dimnames(y)[[1]]
	
	temp <- overlapSets(dataMatrix,tdataMatrix)
	dataMatrix <- temp$x
	tdataMatrix <- temp$y
	sfeatureNames<-row.names(dataMatrix)
	
	# standardize both sets
	if(std){
		dataMatrix<-standardize(dataMatrix)
		tdataMatrix<-standardize(tdataMatrix)
	}
	
	thisClass <- as.vector(classes)
	nclasses<-nlevels(as.factor(thisClass))
	classLevels<-levels(as.factor(thisClass))
	for(j in 1:nclasses){
		thisClass[thisClass==classLevels[j]] <- j
	}
	thisClass<-as.numeric(thisClass)
	dm <- dataMatrix[,!(is.na(thisClass))]
	thisClass <- thisClass[!(is.na(thisClass))]
		
	activeGenes <- min(activeGenes,as.numeric(sprintf("%.0f",dim(dataMatrix)[1]/nclasses)))
	
	trained<-trainClanc(dm, thisClass, sfeatureNames)
	trainedModel<-buildClanc(dm, thisClass, levels(as.factor(thisClass)), trained, activeGenes)
	prediction<-predictClanc(tdataMatrix, sfeatureNames, trainedModel)
	
	testData<-tdataMatrix[row.names(tdataMatrix) %in% trainedModel$geneNames,]	
		
	for(j in 1:nclasses){
		prediction[prediction==j] <- classLevels[j]
	}
		
	prediction<-cbind(tsampleNames,prediction)
	
	temp <- overlapSets(trainedModel$cntrds,tdataMatrix)
	centroids <- temp$x
	testData <- temp$y
	sfeatureNames<-row.names(trainedModel$cntrds)

	distMat<-matrix(ncol=nclasses,nrow=dim(testData)[2])
	for(i in 1:dim(testData)[2]){
		distMat[i,]<-distClanc(testData[,i], centroids, trainedModel$pooledSD, trainedModel$prior)
	}
	dimnames(distMat)<-list(tsampleNames,classLevels)

	dimnames(centroids)<-list(dimnames(centroids)[[1]],classLevels)
	return(list(predictions=prediction,testData=testData,distances=distMat,centroids=centroids))
}

pamPredict<-function(x,classes,y,delta,priors="equal",std=F){

	dataMatrix<-x
	features<- dim(x)[2]
	samples<- dim(x)[1]
	sampleNames<- dimnames(x)[[2]]
	featureNames<- dimnames(x)[[1]]
	
	#parse the test file - same as train file but no rows of classes
	tdataMatrix<-y
	tfeatures<- dim(y)[2]
	tsamples<- dim(y)[1]
	tsampleNames<- dimnames(y)[[2]]
	tfeatureNames<- dimnames(y)[[1]]
	
	temp <- overlapSets(dataMatrix,tdataMatrix)
	dataMatrix <- temp$x
	tdataMatrix <- temp$y
	sfeatureNames<-row.names(dataMatrix)
	
	# standardize both sets
	if(std){
		dataMatrix<-standardize(dataMatrix)
		tdataMatrix<-standardize(tdataMatrix)
	}
	
	# transform the supervised variable to an integer
	thisClass <- as.vector(classes)
	nClasses<-nlevels(as.factor(thisClass))
	classLevels<-levels(as.factor(thisClass))
	for(j in 1:nClasses){
		thisClass[thisClass==classLevels[j]] <- j
	}
	thisClass<-as.numeric(thisClass)
	dm <- dataMatrix[,!(is.na(thisClass))]
	thisClass <- thisClass[!(is.na(thisClass))]
		
	if(priors=="class"){
		priors<-NULL
	}
	if(priors=="equal"){
		priors<-rep(1/nClasses,nClasses)
	}
	
	trainData <- list(x=dm,y=thisClass, geneid=featureNames, genenames=featureNames)
	mytrain <- pamr.train(trainData,prior=priors)
	#new.scales <- pamr.adaptthresh(mytrain)
	#mytrain <- pamr.train(trainData, threshold.scale=new.scales)
	prediction<-pamr.predict(mytrain, tdataMatrix, delta, type="class")
	prediction2<-as.numeric(prediction)
	for(j in 1:nClasses){
		prediction2[prediction==j] <- classLevels[j]
	}
	
	# output the genes and shrunken centroids
	genes<-pamr.listgenes(mytrain, trainData,  delta, genenames=FALSE)
	centroids<-mytrain$centroids
	centroids<-centroids[row.names(centroids) %in% genes[,1],]
	dimnames(centroids)[[2]]<-classLevels	
	
	centroids<-centroids[sort.list(dimnames(centroids)[[1]]),]
	testData<-tdataMatrix[row.names(tdataMatrix) %in% dimnames(centroids)[[1]],]	
	testData<-testData[sort.list(row.names(testData)),]

	distances<-matrix(ncol=nClasses,nrow=dim(tdataMatrix)[2])
	correlations<-matrix(ncol=nClasses,nrow=dim(tdataMatrix)[2])
	for(j in 1:nClasses){
		distances[,j]<-dist(t(cbind(centroids[,j],testData)))[1:dim(tdataMatrix)[2]]
		correlations[,j]<-cor(cbind(centroids[,j],testData))[2:(dim(tdataMatrix)[2]+1)]
	}
	dimnames(distances)<-list(tsampleNames,classLevels)
	dimnames(correlations)<-list(tsampleNames,classLevels)

	return(list(predictions=prediction2,testData=testData,distances=distances,centroids=centroids,correlations=correlations))
}


sspCV<-function(x,classes,geneSetSizes=c(seq(1,20,by=1)),folds=10){

	dataMatrix<-x
	features<- dim(dataMatrix)[1]
	samples<- dim(dataMatrix)[2]
	sampleNames<- names(dataMatrix)
	geneNames<-row.names(dataMatrix)
	
	outTables<-list()
	# loop and do CV for each set of class labels

	# format the dependent variable
	tclass <- as.vector(classes)
	nclasses<-nlevels(as.factor(tclass))
	classLevels<-levels(as.factor(tclass))
	for(j in 1:nclasses){
		tclass[tclass==classLevels[j]] <- j
	}
	tclass<-as.numeric(tclass)
	dm <- dataMatrix[,!(is.na(tclass))]
	tclass <- tclass[!(is.na(tclass))]
	
	# run CV
	cvResults<-cvSSP(dm, tclass, cvIdx = balancedFolds(tclass, folds), active = geneSetSizes)
	
	outTable<-cbind(geneSetSizes,geneSetSizes,cvResults$overallError,cvResults$classError)
	dimnames(outTable)<-list(NULL,c("nGenes per class","nGenes","Overall Error",paste(classLevels,"Error")))
		
	Return(outTable)
}

cvSSP<-function(dm, tclass, cvIdx, active){
	
	nClasses<-max(as.numeric(tclass))

	overallAccuracy<-rep(0,length(active))
	folds<-length(cvIdx)
	classAccuracy<-matrix(0,nrow=length(active),ncol=nClasses)
	for(i in 1:folds){	
	
		dmtrain<-dm[,-1*(cvIdx[[i]])]
		dmtest<-dm[,cvIdx[[i]]]
		trainLabels<-tclass[-1*(cvIdx[[i]])]
		testLabels<-tclass[cvIdx[[i]]]
		
		scores<-apply(dmtrain,1,bwss,trainLabels)
		trainscores<-vector()	
		for(j in 1:dim(dmtrain)[1]){			
			trainscores[j]<-scores[[row.names(dmtrain)[j]]]$bss / scores[[row.names(dmtrain)[j]]]$wss
		}
		
		dmtrain<-dmtrain[sort.list(trainscores,decreasing=T),]
		dmtest<-dmtest[sort.list(trainscores,decreasing=T),]		
		
		for(k in 1:length(active)){
			nGenes<-active[k]
	
			dmtrainG<-dmtrain[1:nGenes,]
			dmtestG<-dmtest[1:nGenes,]
	  	
			centroids<-matrix(,nrow=nGenes,ncol=nClasses)
			for(j in 1:nClasses){
				centroids[,j]<-apply(dmtrainG[,trainLabels==j],1,mean)
			}
			
			distances<-matrix(ncol=nClasses,nrow=dim(dmtest)[2])
			for(j in 1:nClasses){
				distances[,j]<-dist(t(cbind(centroids[,j],dmtestG)))[1:dim(dmtestG)[2]]
			}

# could simplify this with the 'match' function
			classes<-rep(1,dim(dmtestG)[2])
			for(j in 2:nClasses){
				classes<-cbind(classes,rep(j,dim(dmtestG)[2]))
			}
			temp<-apply(distances,1,min)
			temp<-distances==temp
			classes<-classes*temp
			overallAccuracy[k]<-overallAccuracy[k]+sum(apply(classes,1,sum)==testLabels)/length(testLabels)/folds

			for(j in 1:nClasses){
				if(length(testLabels[testLabels==j])==1){
					classAccuracy[k,j]<-classAccuracy[k,j]+sum(classes[testLabels==j,]==testLabels[testLabels==j])/length(testLabels[testLabels==j])/folds
				}else{
					classAccuracy[k,j]<-classAccuracy[k,j]+sum(apply(classes[testLabels==j,],1,sum)==testLabels[testLabels==j])/length(testLabels[testLabels==j])/folds
				}
			}
		}
		print(i)		
	}
	
	return(list(overallError=(1-overallAccuracy),classError=(1-classAccuracy)))
	
}

sspPredict<-function(x,classes="",y,nGenes="",priors="equal",std=F,distm="euclidean",centroids=F){

	dataMatrix<-x
	features<- dim(x)[1]
	samples<- dim(x)[2]
	sampleNames<- dimnames(x)[[2]]
	featureNames<- dimnames(x)[[1]]
	
	#parse the test file - same as train file but no rows of classes
	tdataMatrix<-y
	tfeatures<- dim(y)[1]
	tsamples<- dim(y)[2]
	tsampleNames<- dimnames(y)[[2]]
	tfeatureNames<- dimnames(y)[[1]]
	
	#dimnames(tdataMatrix)[[2]]<-paste("x",seq(1,471))
	temp <- overlapSets(dataMatrix,tdataMatrix)
	dataMatrix <- temp$x
	tdataMatrix <- temp$y
	sfeatureNames<-row.names(dataMatrix)
	
	# standardize both sets
	if(std){
		dataMatrix<-standardize(dataMatrix)
		tdataMatrix<-standardize(tdataMatrix)
	}
	
	if(!centroids){
		thisClass <- as.vector(classes[,1])
		nClasses<-nlevels(as.factor(thisClass))
		classLevels<-levels(as.factor(thisClass))
		for(j in 1:nClasses){
			thisClass[thisClass==classLevels[j]] <- j
		}
		thisClass<-as.numeric(thisClass)
		dataMatrix <- dataMatrix[,!(is.na(thisClass))]
		thisClass <- thisClass[!(is.na(thisClass))]
	
		scores<-apply(dataMatrix,1,bwss,thisClass)
		trainscores<-vector()	
		for(j in 1:dim(dataMatrix)[1]){			
			trainscores[j]<-scores[[row.names(dataMatrix)[j]]]$bss / scores[[row.names(dataMatrix)[j]]]$wss
		}
		
		dataMatrix<-dataMatrix[sort.list(trainscores,decreasing=T),]
		tdataMatrix<-tdataMatrix[sort.list(trainscores,decreasing=T),]	
		
		if(nGenes==""){
			nGenes<-dim(dataMatrix)[1]
		}
		print(paste("Number of genes used:",nGenes))
		
		dataMatrix<-dataMatrix[1:nGenes,]
		tdataMatrix<-tdataMatrix[1:nGenes,]

		centroids<-matrix(,nrow=nGenes,ncol=nClasses)
		for(j in 1:nClasses){
			centroids[,j]<-apply(dataMatrix[,thisClass==j],1,mean)
		}
		dimnames(centroids)<-list(row.names(dataMatrix),NULL)
		
	}else{
		nGenes<-dim(dataMatrix)[1]
		print(paste("Number of genes used:",nGenes))
		centroids<-dataMatrix
		nClasses<-dim(centroids)[2]
		classLevels<-dimnames(centroids)[[2]]
	}
	
	distances<-matrix(ncol=nClasses,nrow=dim(tdataMatrix)[2])
	for(j in 1:nClasses){
		if(distm=="euclidean"){
			distances[,j]<-dist(t(cbind(centroids[,j],tdataMatrix)))[1:(dim(tdataMatrix)[2])]
		}
		if(distm=="correlation" | distm=="pearson"){
			distances[,j]<-t(-1*cor(cbind(centroids[,j],tdataMatrix),use="pairwise.complete.obs"))[2:(dim(tdataMatrix)[2]+1)]
		}
		if(distm=="spearman"){
			distances[,j]<-t(-1*cor(cbind(centroids[,j],tdataMatrix),method="spearman",use="pairwise.complete.obs"))[2:(dim(tdataMatrix)[2]+1)]
		}
	}
	
	scores<-apply(distances,1,min)
	prediction<-vector(length=tsamples)
	for(i in 1:tsamples){
		prediction[i]<-classLevels[match(scores[i],distances[i,])]
	}
	names(prediction)<-tsampleNames
	
	return(list(predictions=prediction,testData=tdataMatrix,distances=distances,centroids=centroids))
}


kmPlot<-function(x,event,stime,varName="",ymin=0,lineColors=NA,nclasses=NA,overall=F,ylab="Probability of Event"){
	
	mainLabel <- varName
	event<-as.numeric(as.vector(event))
	stime<-as.numeric(as.vector(stime))
	if(is.numeric(x)){
		tclass <- factor(x[!is.na(x)])
		event <- event[!is.na(x)]
		stime <- stime[!is.na(x)]
	}else{
		tclass <- factor(x[x!=""])
		event <- event[x!=""]
		stime <- stime[x!=""]
	}
	
	if(is.na(nclasses)){
		nclasses<-nlevels(tclass)
	}
	
	if(length(lineColors)<=1){
		lineColors<-seq(1,nclasses)
	}
	
	y<-survfit(Surv(stime, event)~tclass)
	plot(y,col=lineColors,main=mainLabel,ylim=c(ymin,1),ylab=ylab)
	
	if(overall==T){
		yo<-survfit(Surv(stime, event)~rep(1,length(tclass)))
		lines(yo,col=1,lwd=2)
		legend("bottomleft",legend=c(levels(tclass),"All"),col=c(lineColors,1),bty="n",lty=rep(1,(nclasses+1)),lwd=c(rep(1,nclasses),2.5))
	}else{
		legend("bottomleft",legend=levels(tclass),col=lineColors,bty="n",lty=rep(1,nclasses))
	}
	
	pvalue<-1-pchisq(survdiff(Surv(stime, event)~tclass)$chisq,nclasses-1)
	legend("bottomright",legend=paste("Log Rank p=",signif(pvalue,3),sep=""),cex=1.2,,bty="n")
	
}



edge<-function(x,classes=NA,times,permutations=500){

	#########################################################################
	#       This function computes p-values for temporal differential       #
	#         expression.                                                   #
	#       Author(s): John Storey, Jeffrey Leek                            #
	#       Language: R                                                     #
	#       References:                                                     #
	#         Storey, J.D., Leek, J.T., Xiao, W., Dai, J.Y., Davis R.W.,    #
	#         A Significance Method for Time Course Microarray Experiments  #
	#         Applied to Two Human Studies, University of Washington        #
	#         Biotatistics Department Technical report (2004),              #
	#         http://faculty.washington.edu/~jstorey/                       #
	#       Arguments:                                                      #
	#         dat - The m x n matrix of expression values.                  #
	#         tme - The time covariate vector.                              #
	#         grp - If appropriate, the vector of group labels.             #
	#         ind - If appropriate, the vector of individual labels.        #
	#         B - The number of null permutations.                          #
	#         null - The null hypothesis that is to be tested.              #
	#         dfo - The degrees of freedom for the spline.                  #
	#         intercept - Specifies if an intercept is to be included.      #
	#         outfile - If the results should be written to a file (not     #
	#           used in the GUI).                                           #
	#         basis - The spline basis, either natural cubic spline or      #
	#           polynomial spline.                                          #
	#         eps - likelihood based convergence criterion when EM          #
	#           algorithm is necessary                                      #
	#         seed - option to set seed for random bootstrap samples        #
	#                                                                       #
	#########################################################################
	#timex <- function(dat, tme, grp = NULL, ind = NULL, B = 100,
	#                  null = c("curve", "linear", "flat"), dfo = NULL,
	#                  intercept = TRUE, outfile = NULL,
	#                  basis = c("ncs", "poly"), eps = .05,
	#                  seed = NULL, match=NULL,
	#                  updatefunc=NULL) {
	

 	features<- dim(x)[1]
 	samples<- dim(x)[2]
 	sampleNames<- names(x)
 	geneNames<-row.names(x)
 	x<-apply(x,2,as.numeric)
 	tempTable<-matrix()
 	t0normalize<-F
 
 
	tempTable<-matrix()
	if(is.data.frame(classes)){
		rowsOfClasses <- dim(classes)[2]
		classes<-t(classes)
	
		for(i in 1:rowsOfClasses){
		
			gr.labels <- as.numeric(as.vector(t(classes[i,])))
		
			dm <- x[,!(is.na(gr.labels))]
			gr.times <- times[!(is.na(gr.labels))]
			firstTime <- min(as.numeric(as.vector(t(gr.times))))
			lastTime <- max(as.numeric(as.vector(t(gr.times))))
			names(dm) <- sampleNames[!(is.na(gr.labels))]
			gr.labels <- gr.labels[!(is.na(gr.labels))]
		
			c1.lastTime <- apply(dm[,gr.times==lastTime & gr.labels==1],1,mean)
			c1.firstTime <- apply(dm[,gr.times==firstTime & gr.labels==1],1,mean)
			c1.fc <- c1.lastTime - c1.firstTime
			c2.lastTime <- apply(dm[,gr.times==lastTime & gr.labels==2],1,mean)
			c2.firstTime <- apply(dm[,gr.times==firstTime & gr.labels==2],1,mean)
			c2.fc <- c2.lastTime - c2.firstTime
			fc <- c2.fc - c1.fc
		
			if(t0normalize){
				dm[,gr.labels==1] <- dm[,gr.labels==1] - c1.firstTime
				dm[,gr.labels==2] <- dm[,gr.labels==2] - c2.firstTime
			}
			
			if(firstTime > 1)
				gr.times <- gr.times - 1
		
			result <- timex(dm, as.numeric(as.vector(t(gr.times))), grp = gr.labels, ind = NULL, B = permutations,
		                  null = "curve", dfo = NULL,
		                  intercept = FALSE, outfile = NULL,
		                  basis = "ncs", eps = .05,
		                  seed = NULL, match=NULL,
		                  updatefunc=NULL)
		
			fdrs <- qvalue(result$p)$qvalues
		
			temp <- cbind(c1.fc, c2.fc, fc, result$lr, result$p, fdrs)
			dimnames(temp) <- list(geneNames,c(paste(row.names(classes)[i],"G1 T",lastTime,"/T",firstTime,sep=""),
									paste(row.names(classes)[i]," G2 T",lastTime,"/T",firstTime,sep=""),
									paste(row.names(classes)[i],"trend ratio"),
									paste(row.names(classes)[i],"LR"),
									paste(row.names(classes)[i],"p-value"),
									paste(row.names(classes)[i],"q-value")))
			temp <- temp[sort.list(row.names(temp)),]
		
			if(i==1){
				tempTable <- temp
				dimnames(tempTable)[[2]] <- c(paste(row.names(classes)[i]," G1 T",firstTime,"/T",lastTime,sep=""),
									paste(row.names(classes)[i]," G2 T",lastTime,"/T",firstTime,sep=""),
									paste(row.names(classes)[i],"trend ratio"),
									paste(row.names(classes)[i],"LR"),
									paste(row.names(classes)[i],"p-value"),
									paste(row.names(classes)[i],"q-value"))
			}else{
				tempTable <- cbind(tempTable,temp)
				tempDim <- dim(tempTable)[2]
				dimnames(tempTable)[[2]][(tempDim-5):tempDim] <- c(paste(row.names(classes)[i]," G1 T",lastTime,"/T",firstTime,sep=""),
									paste(row.names(classes)[i]," G2 T",lastTime,"/T",firstTime,sep=""),
									paste(row.names(classes)[i],"trend ratio"),
									paste(row.names(classes)[i],"LR"),
									paste(row.names(classes)[i],"p-value"),
									paste(row.names(classes)[i],"q-value"))
			}
		}
	}else{
			gr.labels <- as.numeric(as.vector(classes))
		
			dm <- x[,!(is.na(gr.labels))]
			gr.times <- times[!(is.na(gr.labels))]
			firstTime <- min(as.numeric(as.vector(t(gr.times))))
			lastTime <- max(as.numeric(as.vector(t(gr.times))))
			names(dm) <- sampleNames[!(is.na(gr.labels))]
			gr.labels <- gr.labels[!(is.na(gr.labels))]
		
			c1.lastTime <- apply(dm[,gr.times==lastTime & gr.labels==1],1,mean)
			c1.firstTime <- apply(dm[,gr.times==firstTime & gr.labels==1],1,mean)
			c1.fc <- c1.lastTime - c1.firstTime
			c2.lastTime <- apply(dm[,gr.times==lastTime & gr.labels==2],1,mean)
			c2.firstTime <- apply(dm[,gr.times==firstTime & gr.labels==2],1,mean)
			c2.fc <- c2.lastTime - c2.firstTime
			fc <- c2.fc - c1.fc
		
			if(t0normalize){
				dm[,gr.labels==1] <- dm[,gr.labels==1] - c1.firstTime
				dm[,gr.labels==2] <- dm[,gr.labels==2] - c2.firstTime
			}
			
			if(firstTime > 1)
				gr.times <- gr.times - 1
		
			result <- timex(dm, as.numeric(as.vector(t(gr.times))), grp = gr.labels, ind = NULL, B = permutations,
		                  null = "curve", dfo = NULL,
		                  intercept = FALSE, outfile = NULL,
		                  basis = "ncs", eps = .05,
		                  seed = NULL, match=NULL,
		                  updatefunc=NULL)
		
			fdrs <- qvalue(result$p)$qvalues
		
			temp <- cbind(c1.fc, c2.fc, fc, result$lr, result$p, fdrs)
			dimnames(temp) <- list(geneNames,c(paste(row.names(classes)[i],"G1 T",lastTime,"/T",firstTime,sep=""),
									paste(row.names(classes)[i]," G2 T",lastTime,"/T",firstTime,sep=""),
									paste(row.names(classes)[i],"trend ratio"),
									paste(row.names(classes)[i],"LR"),
									paste(row.names(classes)[i],"p-value"),
									paste(row.names(classes)[i],"q-value")))
			temp <- temp[sort.list(row.names(temp)),]
			tempTable <- temp
			dimnames(tempTable)[[2]] <- c(paste(row.names(classes)[i]," G1 T",firstTime,"/T",lastTime,sep=""),
									paste(row.names(classes)[i]," G2 T",lastTime,"/T",firstTime,sep=""),
									paste(row.names(classes)[i],"trend ratio"),
									paste(row.names(classes)[i],"LR"),
									paste(row.names(classes)[i],"p-value"),
									paste(row.names(classes)[i],"q-value"))

	}
	return(tempTable)
}


edgeOneClass<-function(x,classes=NA,times,permutations=500){

 	features<- dim(x)[1]
 	samples<- dim(x)[2]
 	sampleNames<- dimnames(x)[[2]]
 	geneNames<- dimnames(x)[[1]]
 	x<-apply(x,2,as.numeric)
 	rowsOfClasses <- dim(classes)[1]
 	tempTable<-matrix()
 	t0normalize<-F
 	
	for(i in 1:rowsOfClasses){

		gr.labels <- as.vector(as.numeric(t(classes[i,])))
	
		dm <- x[,!(is.na(gr.labels))]
		gr.times <- as.vector(as.numeric(t(times[!(is.na(gr.labels))])))
		firstTime <- min(gr.times) 
		lastTime <- max(gr.times)
		names(dm) <- sampleNames[!(is.na(gr.labels))]
		gr.labels <- gr.labels[!(is.na(gr.labels))]
	
		if(sum(gr.times==lastTime & gr.labels==1)>1){
			c1.lastTime <- apply(dm[,gr.times==lastTime & gr.labels==1],1,mean)
		}else{
			c1.lastTime <- dm[,gr.times==lastTime & gr.labels==1]
		}
		if(sum(gr.times==firstTime & gr.labels==1)>1){
			c1.firstTime <- apply(dm[,gr.times==firstTime & gr.labels==1],1,mean)
		}else{
			c1.firstTime <- dm[,gr.times==firstTime & gr.labels==1]
		}
		c1.fc <- c1.lastTime - c1.firstTime

		
		if(t0normalize){
			dm[,gr.labels==1] <- dm[,gr.labels==1] - c1.firstTime
		}
		
		if(firstTime > 1)
			gr.times <- gr.times - 1
	
		result <- timex(dm, gr.times, ind = NULL, B = permutations,
	                  null = "flat", dfo = NULL,
	                  intercept = FALSE, outfile = NULL,
	                  basis = "ncs", eps = .05,
	                  seed = NULL, match=NULL,
	                  updatefunc=NULL)
	
		fdrs <- qvalue(result$p)$qvalues
	
		temp <- cbind(c1.fc, result$lr, result$p, fdrs)
		dimnames(temp) <- list(geneNames,c(
								paste(row.names(classes)[i],"WT T",lastTime,"/T",firstTime,sep=""),
								paste(row.names(classes)[i],"LR"),
								paste(row.names(classes)[i],"p-value"),
								paste(row.names(classes)[i],"q-value")))
		temp <- temp[sort.list(row.names(temp)),]
	
		if(i==1){
			tempTable <- temp
			dimnames(tempTable)[[2]] <- c(paste(row.names(classes)[i]," WT T",firstTime,"/T",lastTime,sep=""),
								paste(row.names(classes)[i],"LR"),
								paste(row.names(classes)[i],"p-value"),
								paste(row.names(classes)[i],"q-value"))
		}else{
			tempTable <- cbind(tempTable,temp)
			tempDim <- dim(tempTable)[2]
			dimnames(tempTable)[[2]][(tempDim-3):tempDim] <- c(paste(row.names(classes)[i]," WT T",lastTime,"/T",firstTime,sep=""),
								paste(row.names(classes)[i],"LR"),
								paste(row.names(classes)[i],"p-value"),
								paste(row.names(classes)[i],"q-value"))
		}
	}
	return(tempTable)
}


collapseIDs<-function(x,allids=row.names(x),method="mean"){

	allids<-as.vector(allids)
	ids<- levels(as.factor(allids))
	x.col<- matrix(nrow=length(ids), ncol=dim(x)[2])

	if(length(ids)==dim(x)[1]){ 
			dimnames(x)[[1]]<-allids
			return(x) 
	}
	
	for(i in 1:length(ids)){
		if(sum(allids==ids[i])>1){
			indices <- allids==ids[i] 
			if(method=="mean"){
				vals<-apply(x[indices,],2,mean,na.rm=T)
			}
			if(method=="sum"){ #JYL 12/12/13
				vals<-apply(x[indices,],2,sum,na.rm=T)#JYL 12/12/13
			}#JYL 12/12/13
			if(method=="median"){
				vals<-apply(x[indices,],2,median,na.rm=T)
			}
			if(method=="stdev"){   
				temp<- x[indices,]
				stdevs<- apply(temp,1,sd,na.rm=T)
				vals<- temp[match(max(stdevs),stdevs),]
			}
			if(method=="iqr"){   
				temp<- x[indices,]
				iqrs<- apply(temp,1,function(x){quantile(x,.75,na.rm=T)-quantile(x,.25,na.rm=T)})
				vals<- temp[match(max(iqrs),iqrs),]
			}
			x.col[i,] <- vals
		}else{
			x.col[i,] <- x[allids==ids[i],]
		}
	}

	dimnames(x.col)<- list(ids,dimnames(x)[[2]])
	return(x.col)
	
}


knnCV<-function(x,classes,geneSetSizes=c(seq(1,20,by=1)),nk=3){

	folds <- 10
	dataMatrix<-x
	features<- dim(dataMatrix)[1]
	samples<- dim(dataMatrix)[2]
	sampleNames<- names(dataMatrix)
	geneNames<-row.names(dataMatrix)
	if (geneSetSizes*5 > features || length(geneSetSizes)*5 > features)
	{
         geneSetSizes = floor(features/5)
	}
	outTables<-list()
	# loop and do CV for each set of class labels

	# format the dependent variable
	tclass <- as.vector(t(classes))
	nclasses<-nlevels(as.factor(tclass))
	classLevels<-levels(as.factor(tclass))
	for(j in 1:nclasses){
		tclass[tclass==classLevels[j]] <- j
	}
	tclass<-as.numeric(tclass)
	dm <- dataMatrix[,!(is.na(tclass))]
	tclass <- tclass[!(is.na(tclass))]
	
	# run CV
	cvResults<-knnCV.internal(dm, tclass, cvIdx = balancedFolds(tclass, folds), active = geneSetSizes, nk)
	
	outTable<-cbind(geneSetSizes,geneSetSizes,cvResults$overallError,cvResults$classError)
	dimnames(outTable)<-list(NULL,c("nGenes per class","nGenes","Overall Error",paste(classLevels,"Error")))

	return(outTable)
}

knnPredict<-function(x, classes, y, ngenes, nk){

	thisClass <- as.vector(classes)
	nClasses<-nlevels(as.factor(thisClass))
	classLevels<-levels(as.factor(thisClass))
	for(j in 1:nClasses){
		thisClass[thisClass==classLevels[j]] <- j
	}
	thisClass<-as.numeric(thisClass)
	x <- x[,!(is.na(thisClass))]
	thisClass <- thisClass[!(is.na(thisClass))]
	
	temp<-overlapSets(x,y)
	x<-temp$x
	y<-temp$y
	scores<-apply(x,1,bwss,thisClass)
	trainscores<-vector()	
	for(j in 1:dim(x)[1]){			
		trainscores[j]<-scores[[row.names(x)[j]]]$bss / scores[[row.names(x)[j]]]$wss
	}
	x<-x[sort.list(trainscores,decreasing=T),]
	y<-y[sort.list(trainscores,decreasing=T),]
	x<-x[1:ngenes,]
	y<-y[1:ngenes,]
	
	kpreds<-knn(t(x), t(y), thisClass, k = nk)
	
	return(list(predictions=kpreds, testData = y))
}

knnCV.internal<-function(dm, tclass, cvIdx, active, nk){
	
	nClasses<-max(as.numeric(tclass))

	overallAccuracy<-rep(0,length(active))
	folds<-length(cvIdx)
	classAccuracy<-matrix(0,nrow=length(active),ncol=nClasses)
	for(i in 1:folds){	
	
		dmtrain<-dm[,-1*(cvIdx[[i]])]
		dmtest<-dm[,cvIdx[[i]]]
		trainLabels<-tclass[-1*(cvIdx[[i]])]
		testLabels<-tclass[cvIdx[[i]]]
		
		scores<-apply(dmtrain,1,bwss,trainLabels)
		trainscores<-vector()	
		for(j in 1:dim(dmtrain)[1]){			
			trainscores[j]<-scores[[row.names(dmtrain)[j]]]$bss / scores[[row.names(dmtrain)[j]]]$wss
		}
		
		dmtrain<-dmtrain[sort.list(trainscores,decreasing=T),]
		dmtest<-dmtest[sort.list(trainscores,decreasing=T),]		
		
		for(k in 1:length(active)){
			nGenes<-active[k]
	
			dmtrainG<-dmtrain[1:nGenes,]
			dmtestG<-dmtest[1:nGenes,]
	  	classes<-knn(t(dmtrainG),t(dmtestG),trainLabels, nk)
	  	classes<-as.numeric(classes)
			overallAccuracy[k]<-overallAccuracy[k]+sum(classes==testLabels)/length(testLabels)/folds

			for(j in 1:nClasses){
					classAccuracy[k,j]<-classAccuracy[k,j]+sum(classes[testLabels==j]==testLabels[testLabels==j])/length(testLabels[testLabels==j])/folds
			}
		}
		print(i)
	}
	
	return(list(overallError=(1-overallAccuracy),classError=(1-classAccuracy)))
	
}



pcaEA<-function(x,classes,size=1,showLegend=T,legendloc="topright",mainStr="",startPC=1,stopPC=2,showNames=T,showClasses=F,axisExpansion=0,groupColors=NA){
	
	features<- dim(x)[1]
	samples<- dim(x)[2]
	sampleNames<- dimnames(x)[[2]]
	featureNames<-dimnames(x)[[1]]
	x<-apply(x,2,as.numeric)
	
	#principal components plots
	data.pca<-prcomp(as.matrix(x))

	# Proportion of total variance distributed over 10 first components:
	tmp<-data.pca$sdev[1:10]^2/sum(data.pca$sdev^2)

	gr.labels<-as.vector(t(classes))
	gr.labels.fac<-factor(as.vector(t(classes)),exclude="")
	nlabels<-nlevels(gr.labels.fac)
	legendLabels<-vector()
	legendColors<-vector()
	for(k in 1:nlabels){
		group<-levels(gr.labels.fac)[k]
		legendLabels[k]<-group
		if(length(groupColors)>1){
				gr.labels[gr.labels.fac==group]<-groupColors[k]
				legendColors[k]<-groupColors[k]
		}else{
				gr.labels[gr.labels.fac==group]<-k
				legendColors[k]<-k
		}
	}
	
	if(length(groupColors)==1){
		gr.labels<-as.numeric(gr.labels)
	}
	
	#plot 2pcs by each other
	i<-startPC
	j<-stopPC
	
	#graphing parameters
	par(lab=c(3,4,3))
 	par(mgp=c(.3,.5,.0))
	par(mai=c(.5,.5,.5,.5))
	par(xaxt="n",yaxt="n")

	strM<-mainStr
	strX<-paste("PC",i,paste("(",round(tmp[i],4)*100,"%)",sep=""),sep=" ")
	strY<-paste("PC",j,paste("(",round(tmp[j],4)*100,"%)",sep=""),sep=" ")
	xmin<-min(data.pca$rotation[,i])-abs(axisExpansion*min(data.pca$rotation[,i]))
	xmax<-max(data.pca$rotation[,i])+abs(axisExpansion*max(data.pca$rotation[,i]))
	ymin<-min(data.pca$rotation[,j])-abs(axisExpansion*min(data.pca$rotation[,j]))
	ymax<-max(data.pca$rotation[,j])+abs(axisExpansion*max(data.pca$rotation[,j]))
	plot(data.pca$rotation[,i],data.pca$rotation[,j], xlab=strX, ylab=strY, 
							main=strM, col=gr.labels,cex=size,pch="",xlim=c(xmin,xmax),ylim=c(ymin,ymax))
	if(showNames){
		text(data.pca$rotation[,i],data.pca$rotation[,j],labels=names(data.pca$rotation[,i]),cex=size*.6)
	}else{
		if(showClasses){
				text(data.pca$rotation[,i],data.pca$rotation[,j],labels=gr.labels.fac,cex=size*.6)
		}else{
				points(data.pca$rotation[,i],data.pca$rotation[,j],col=gr.labels,cex=size*1.5,pch=19)
		}
		if(showLegend){
			legend(legendloc,legend=legendLabels,col=legendColors,pch=19,x.intersp=.3,yjust=.5,bty="n",cex=size)
		}
	}
}

quantileScrunch<-function(x){
#scrunch to 1-2 instead of 0-1
	temp<-dimnames(x)
	scrunch <- x<=1
	mins<-min(x[scrunch])
	maxs<-max(x[scrunch])
	x[scrunch]<-((x[scrunch]-mins)/(maxs-mins))+1
	x<-normalize.quantiles(as.matrix(x))
	dimnames(x)<-temp
	return(x)
}



svmCV<-function(x,classes,geneSetSizes=seq(20,1,by= -1),kernel="linear"){
	folds <- 10
	dataMatrix<-x
	features<- dim(dataMatrix)[1]
	samples<- dim(dataMatrix)[2]
	sampleNames<- names(dataMatrix)
	geneNames<-row.names(dataMatrix)

	# format the dependent variable
	tclass <- as.vector(t(classes))
	nclasses<-nlevels(as.factor(tclass))
	classLevels<-levels(as.factor(tclass))
	for(j in 1:nclasses){
		if(j==1){
			tclass[tclass==classLevels[j]] <- -1
		}else{
			if(j==2){
				tclass[tclass==classLevels[j]] <- 1
			}else{
				tclass[tclass==classLevels[j]] <- NA
			}
		}
	}
	tclass<-as.numeric(tclass)
	dm <- t(dataMatrix[,!(is.na(tclass))])
	tclass <- tclass[!(is.na(tclass))]
	
	# run CV
	cvResults<-RSVM(dm, tclass, geneSetSizes, CVtype=10, kernfunc=kernel)
	
	outTable<-cbind(cvResults$ladder,cvResults$ladder,cvResults$Error)
	dimnames(outTable)<-list(NULL,c("nGenes per class","nGenes","Overall Error"))
	rownames(cvResults$SelFreq)<-rownames(x)
	
	return(list(outTable,cvResults$SelFreq))
}

svmPredict<-function(x, classes, y, cvFreq, ngenes, kernel="linear"){

	# format the dependent variable
	tclass <- as.vector(classes)
	nclasses<-nlevels(as.factor(tclass))
	classLevels<-levels(as.factor(tclass))
	for(j in 1:nclasses){
		if(j==1){
			tclass[tclass==classLevels[j]] <- -1
		}else{
			if(j==2){
				tclass[tclass==classLevels[j]] <- 1
			}else{
				tclass[tclass==classLevels[j]] <- NA
			}
		}
	}
	tclass<-as.numeric(tclass)
	dm <- x[,!(is.na(tclass))]
	tclass <- tclass[!(is.na(tclass))]
	
	temp<-overlapSets(dm,y)
	x<-t(temp$x)
	y<-t(temp$y)
	
	gcol<-paste("Lev_",ngenes,sep="")
	vars<-cvFreq[[2]][,gcol]
	
	x<-x[,sort.list(vars,decreasing=T)[1:ngenes]]
	y<-y[,sort.list(vars,decreasing=T)[1:ngenes]]

	svmres <- svm(x, factor(tclass), scale=F, type="C-classification", kernel=kernel )
	svmpred <- predict(svmres, y )
	
	return(list(predictions=svmpred, testData = t(y)))
}

### R-code for R-SVM
### use leave-one-out / Nfold or bootstrape to permute data for external CV
### build SVM model and use mean-balanced weight to sort genes on training set
### and recursive elimination of least important genes
### author: Dr. Xin Lu, Research Scientist
###   Biostatistics Department, Harvard School of Public Health


## R-SVM core code
## input:
##    x: row matrix of data
##    y: class label: 1 / -1 for 2 classes
##    CVtype: 
##        integer: N fold CV
##        "LOO":    leave-one-out CV
##        "bootstrape": bootstrape CV
##    CVnum:   number of CVs
##        LOO: defined as sample size
##        Nfold and bootstrape:  user defined, default as sample size
## output: a named list
##    Error: a vector of CV error on each level
##    SelFreq: a matrix for the frequency of each gene being selected in each level
##             with each column corresponds to a level of selection
##             and each row for a gene
##          The top important gene in each level are those high-freqent ones
RSVM <- function(x, y, ladder, CVtype, CVnum=0, kernfunc="linear")
{
    ## check if y is binary response
    Ytype <- names(table(y))
    if( length(Ytype) != 2) 
    {
        print("ERROR!! RSVM can only deal with 2-class problem")
        return(0)
    }

    ## class mean
    m1 <- apply(x[ which(y==Ytype[1]), ], 2, mean)
    m2 <- apply(x[ which(y==Ytype[2]), ], 2, mean)
    md <- m1-m2

    yy <- vector( length=length(y))
    yy[which(y==Ytype[1])] <- 1
    yy[which(y==Ytype[2])] <- -1        
    y <- yy

    ## check ladder
    if( min(diff(ladder)) >= 0 )
    {
        print("ERROR!! ladder must be monotonously decreasing")
        return(0);
    }
    
    if( ladder[1] != ncol(x) )
    {
        ladder <- c(ncol(x), ladder)
    }

    nSample <- nrow(x)
    nGene   <- ncol(x)
    SampInd <- seq(1, nSample)

    if( CVtype == "LOO" )
    {
        CVnum <- nSample
    } else
    {
        if( CVnum == 0 )
        {
            CVnum <- nSample
        }
    }
    
    ## vector for test error and number of tests
    ErrVec <- vector( length=length(ladder))
    names(ErrVec) <- paste("Lev_", ladder, sep="")
    nTests <- 0

    SelFreq <- matrix( 0, nrow=nGene, ncol=length(ladder))
    colnames(SelFreq) <- paste("Lev_", ladder, sep="")
    
    ## for each CV    
    for( i in 1:CVnum )
    {
    
        ## split data
        if( CVtype == "LOO" )
        {
            TestInd <- i
            TrainInd <- SampInd[ -TestInd]
        } else
        {
            if( CVtype == "bootstrape" )
            {
                TrainInd <- sample(SampInd, nSample, replace=T )
                TestInd <- SampInd[ which(!(SampInd %in% TrainInd ))]
            } else
            {
                ## Nfold
                TrainInd <- sample(SampInd, nSample*(CVtype-1)/CVtype )
                TestInd <- SampInd[ which(!(SampInd %in% TrainInd ))]
            }
        }

        nTests <- nTests + length(TestInd)
        
        ## in each level, train a SVM model and record test error
        xTrain <- x[TrainInd, ]
        yTrain <- y[TrainInd]
       
        xTest  <- x[TestInd,]
        yTest  <- y[TestInd]

        ## index of the genes used in the 
        SelInd <- seq(1, nGene)
        for( gLevel in 1:length(ladder) )
        {
            ## record the genes selected in this ladder
            SelFreq[SelInd, gLevel] <- SelFreq[SelInd, gLevel] +1
            
            ## train SVM model and test error
             svmres <- svm(xTrain[, SelInd], yTrain, scale=F, type="C-classification", kernel=kernfunc )
             if( CVtype == "LOO" )
             {
                 svmpred <- predict(svmres, matrix(xTest[SelInd], nrow=1) )
             } else
             {
                 svmpred <- predict(svmres, xTest[, SelInd] )
             }
             ErrVec[gLevel] <- ErrVec[gLevel] + sum(svmpred != yTest )
             
            ## weight vector
             W <- t(svmres$coefs*yTrain[svmres$index]) %*% svmres$SV * md[SelInd]
             rkW <- rank(W)
             
             if( gLevel < length(ladder) )
             {
                SelInd <- SelInd[which(rkW > (ladder[gLevel] - ladder[gLevel+1]))]
             }
        }

    }

    ret <- list(ladder=ladder, Error=ErrVec/nTests, SelFreq=SelFreq)
    
}


dldaCV<-function(x,classes,geneSetSizes=c(seq(2,20,by=1))){

	folds <- 10
	dataMatrix<-x
	features<- dim(dataMatrix)[1]
	samples<- dim(dataMatrix)[2]
	sampleNames<- names(dataMatrix)
	geneNames<-row.names(dataMatrix)
	
	# format the dependent variable
	tclass <- as.vector(classes)
	nclasses<-nlevels(as.factor(tclass))
	classLevels<-levels(as.factor(tclass))
	for(j in 1:nclasses){
		tclass[tclass==classLevels[j]] <- j
	}
	tclass<-as.numeric(tclass)
	dm <- dataMatrix[,!(is.na(tclass))]
	tclass <- tclass[!(is.na(tclass))]
	
	# run CV
	cvResults<-dldaCV.internal(dm, tclass, cvIdx = balancedFolds(tclass, folds), active = geneSetSizes)
	
	outTable<-cbind(geneSetSizes,geneSetSizes,cvResults$overallError,cvResults$classError)
	dimnames(outTable)<-list(NULL,c("nGenes per class","nGenes","Overall Error",paste(classLevels,"Error")))

	return(outTable)
}

dldaPredict<-function(x, classes, y, ngenes){

	thisClass <- as.vector(t(classes[1,]))
	nClasses<-nlevels(as.factor(thisClass))
	classLevels<-levels(as.factor(thisClass))
	for(j in 1:nClasses){
		thisClass[thisClass==classLevels[j]] <- j
	}
	thisClass<-as.numeric(thisClass)
	x <- x[,!(is.na(thisClass))]
	thisClass <- thisClass[!(is.na(thisClass))]
	
	temp<-overlapSets(x,y)
	x<-temp$x
	y<-temp$y
	scores<-apply(x,1,bwss,thisClass)
	trainscores<-vector()	
	for(j in 1:dim(x)[1]){			
		trainscores[j]<-scores[[row.names(x)[j]]]$bss / scores[[row.names(x)[j]]]$wss
	}
	x<-x[sort.list(trainscores,decreasing=T),]
	y<-y[sort.list(trainscores,decreasing=T),]
	x<-x[1:ngenes,]
	y<-y[1:ngenes,]
	
	kpreds<-dlda(t(x), t(y), thisClass)
	
	return(list(predictions=kpreds, testData = y))
}

dwdCV<-function(x,y,cvFolds=10,nGenes=50){
	
	nfeatures<-dim(x)[1]
	x<-x[,!is.na(as.vector(t(y)))]
	y<-as.vector(t(y[!is.na(y)]))
	nsamples<-dim(x)[2]
	
	nClasses<-nlevels(as.factor(y))
	if(nClasses>2){
		classLevels<-levels(as.factor(y))
		classMat<-matrix(nrow=nClasses,ncol=nsamples)
		dimnames(classMat)<-list(classLevels,dimnames(x)[[2]])
		for(j in 1:nClasses){
			classMat[j,y==classLevels[j]] <- 1
			classMat[j,y!=classLevels[j]] <- 2
		}
	}else{
		classLevels<-levels(as.factor(y))
		classMat<-matrix(nrow=1,ncol=nsamples)
		dimnames(classMat)<-list(NA,dimnames(x)[[2]])
		classMat[1,y==classLevels[1]] <- 1
		classMat[1,y!=classLevels[1]] <- 2
	}	
	
	#now CV
	cvRuns<-dim(classMat)[1]
	classPreds<-matrix(nrow=nsamples,ncol=cvRuns)
	folds<-balancedFolds(y,cvFolds)
	for(i in 1:cvRuns){
		dwdout<-dwdCV.inner(x,classMat[i,],folds,nGenes)
		classPreds[,i]<-dwdout$cvPreds
	}
	
	scores<-apply(classPreds,1,max)
	prediction<-vector(length=nsamples)
	for(i in 1:nsamples){
		prediction[i]<-classLevels[match(scores[i],classPreds[i,])]
	}
	
	overallAccuracy=sum(prediction==y)/nsamples
	
	classAccuracy<-matrix(0,nrow=1,ncol=nClasses)
	dimnames(classAccuracy)[[2]]<-classLevels
	for(j in 1:nClasses){
		classAccuracy[1,j]<-sum(prediction[y==classLevels[j]]==y[y==classLevels[j]])/sum(y==classLevels[j])
	}
	
	return(list(overallError=(1-overallAccuracy),classError=(1-classAccuracy)))
}
		

dwdCV.inner<-function(x,y,folds,nGenes){

	nfeatures<-dim(x)[1]
	fullDataMatrix<-x
	classes<-y
	nsamples<-dim(fullDataMatrix)[2]
	
	cvPreds<-vector(length=nsamples)

	geneWeights<-matrix(nrow=nfeatures, ncol=length(folds))
	dimnames(geneWeights)<-list(row.names(x),NULL)

	for(cv in 1:length(folds)){
		testDataMatrix<-as.matrix(fullDataMatrix[,folds[[cv]]])
		dataMatrix<-fullDataMatrix[,(-1*folds[[cv]])]
		cvClasses<-as.numeric(classes[(-1*folds[[cv]])])
		
		# create a test set
		# need to implement CV here
		dwddir<-dwd(dataMatrix,cvClasses)
		geneWeights[,cv]<-dwddir
		
		# reduce the gene set and re-compute DWD direction
		dataMatrix<-dataMatrix[row.names(dataMatrix) %in% names(sort(abs(geneWeights[,cv]),decreasing=T))[1:nGenes],]
		dwddir<-dwd(dataMatrix,cvClasses)
		testDataMatrix<-testDataMatrix[row.names(testDataMatrix) %in% row.names(dataMatrix),]
		
		# project samples in the dwd direction
		z <- t(dataMatrix) %*% as.vector(dwddir)
		zz <- t(testDataMatrix) %*% as.vector(dwddir)
		
		# now logistic regression
		
		# format data
		trainingData <- z
		testData <- zz
		classes2<-cvClasses
		classes2[cvClasses==2]<-0
		trainingSet <- list(x=trainingData,y=classes2)
		testSet <- list(x=testData,y=rep(0,dim(testDataMatrix)[2]))
		
		# run logistic regression
		lr.out <- glm(y~x, family=binomial(logit), data=trainingSet)

		cvPreds[folds[[cv]]]<-predict(lr.out,testSet,type="response")
		print(cv)
	}	
	names(cvPreds)<-dimnames(fullDataMatrix)[[2]]
	return(list(cvPreds=cvPreds,weights=geneWeights,classOrder=classLevels))
}


dwdPredict<-function(x,classes,y,std=F){
	
	nfeatures<-dim(x)[1]
	x<-x[,!is.na(as.vector(t(classes)))]
	classes<-as.vector(t(classes[!is.na(classes)]))
	nsamples<-dim(x)[2]
	ntestSamples<-dim(y)[2]
	
	#dimnames(tdataMatrix)[[2]]<-paste("x",seq(1,471))
	temp <- overlapSets(x,y)
	x <- temp$x
	y <- temp$y
	
	# standardize both sets
	if(std){
		x<-standardize(x)
		y<-standardize(y)
	}
	
	nClasses<-nlevels(as.factor(classes))
	if(nClasses>2){
		classLevels<-levels(as.factor(classes))
		classMat<-matrix(nrow=nClasses,ncol=nsamples)
		dimnames(classMat)<-list(classLevels,dimnames(x)[[2]])
		for(j in 1:nClasses){
			classMat[j,classes==classLevels[j]] <- 1
			classMat[j,classes!=classLevels[j]] <- 2
		}
	}else{
		classLevels<-levels(as.factor(classes))
		classMat<-matrix(nrow=1,ncol=nsamples)
		#dimnames(classMat)<-list(classLevels,dimnames(x)[[2]])
		classMat[1,classes==classLevels[1]] <- 1
		classMat[1,classes!=classLevels[1]] <- 2
		nClasses<-1
	}	
	
	#now prediction for each class
	classPreds<-matrix(nrow=ntestSamples,ncol=nClasses)
	classProbs<-matrix(nrow=ntestSamples,ncol=nClasses)
	for(i in 1:nClasses){
		dwddir<-dwd(x,classMat[i,])
		geneWeights<-dwddir
			
		# project samples in the dwd direction
		z <- t(x) %*% as.vector(dwddir)
		zz <- t(y) %*% as.vector(dwddir)
		
		# now logistic regression
		
		# format data
		trainingData <- z
		testData <- zz
		classes2<-classMat[i,]
		classes2[classMat[i,]==2]<-0
		trainingSet <- list(x=trainingData,y=classes2)
		testSet <- list(x=testData,y=rep(0,dim(y)[2]))
		
		# run logistic regression
		lr.out <- glm(y~x, family=binomial(logit), data=trainingSet)

		classPreds[,i]<-predict(lr.out,testSet)
		classProbs[,i]<-predict(lr.out,testSet,type="response")
	}
	
	if(nClasses>1){
		scores<-apply(classPreds,1,max)
		prediction<-vector(length=ntestSamples)
		for(i in 1:ntestSamples){
			prediction[i]<-classLevels[match(scores[i],classPreds[i,])]
		}
		dimnames(classPreds)<-list(dimnames(y)[[2]],classLevels)
	}else{
		prediction<-vector(length=ntestSamples)
		prediction[classPreds[,1]>0.5]<-classLevels[1]
		prediction[classPreds[,1]<0.5]<-classLevels[2]
	}
	
	names(prediction)<-dimnames(y)[[2]]
	return(list(predictions=prediction,probability=classProbs,linear=classPreds))
}



dldaCV.internal<-function(dm, tclass, cvIdx, active){
	
	nClasses<-max(as.numeric(tclass))

	overallAccuracy<-rep(0,length(active))
	folds<-length(cvIdx)
	classAccuracy<-matrix(0,nrow=length(active),ncol=nClasses)
	for(i in 1:folds){	
	
		dmtrain<-dm[,-1*(cvIdx[[i]])]
		dmtest<-dm[,cvIdx[[i]]]
		trainLabels<-tclass[-1*(cvIdx[[i]])]
		testLabels<-tclass[cvIdx[[i]]]
		
		scores<-apply(dmtrain,1,bwss,trainLabels)
		trainscores<-vector()	
		for(j in 1:dim(dmtrain)[1]){			
			trainscores[j]<-scores[[row.names(dmtrain)[j]]]$bss / scores[[row.names(dmtrain)[j]]]$wss
		}
		
		dmtrain<-dmtrain[sort.list(trainscores,decreasing=T),]
		dmtest<-dmtest[sort.list(trainscores,decreasing=T),]		
		
		for(k in 1:length(active)){
			nGenes<-active[k]
	
			dmtrainG<-dmtrain[1:nGenes,]
			dmtestG<-dmtest[1:nGenes,]
	  	classes<-dlda(t(dmtrainG),t(dmtestG),trainLabels)
	  	classes<-as.numeric(classes)
			overallAccuracy[k]<-overallAccuracy[k]+sum(classes==testLabels)/length(testLabels)/folds

			for(j in 1:nClasses){
					classAccuracy[k,j]<-classAccuracy[k,j]+sum(classes[testLabels==j]==testLabels[testLabels==j])/length(testLabels[testLabels==j])/folds
			}
		}
		print(i)
	}
	
	return(list(overallError=(1-overallAccuracy),classError=(1-classAccuracy)))
	
}

dwd<-function(x,y){

# requires a very specific directory structure in the root directory on C:
# file names and other run time requirements are necessitated by the exe

	y<-as.numeric(as.vector(t(y)))
	x<-rbind(y,x)
	
	cwd<-getwd()
	setwd("C:/DWD/DWDdata")
	write.table(x,"DWD_Input.txt",sep="\t",col.names=F,row.names=F)

	twd<-getwd()
	write(twd,"fileDir.txt")

	y<-0
	write(y,"DWD_MeanAdjustType.txt")

	setwd("C:/DWD/lib")
	system("BatchAdjustSM.exe",invisible=T)

	setwd("C:/DWD/DWDdata")
	unlink("fileDir.txt")
	unlink("DWD_MeanAdjustType.txt")
	unlink("DWD_Input.txt")

	setwd("C:/DWD/DWDdata")
	dwdVec<-read.table("DWD_Vec.txt",sep="\t",header=F)

	setwd(cwd)

	return(t(dwdVec))
}


dwdadjust<-function(x,y,dwddir){

	x1<-x[,y==1]
	x2<-x[,y==2]

	traindwdp <- t(x1) %*% as.vector(dwddir)
	traindwdpm <- mean(traindwdp) * dwddir
	testdwdp <- t(x2) %*% as.vector(dwddir)
	testdwdpm <- mean(testdwdp) * dwddir

	trainadjmat <- matrix(nrow=dim(x1)[1],ncol=dim(x1)[2])
	for(i in 1:dim(x1)[2]){
			trainadjmat[,i]<-x1[,i]-t(traindwdpm)
	}
	testadjmat <- matrix(nrow=dim(x1)[1],ncol=dim(x2)[2])
	for(i in 1:dim(x2)[2]){
			testadjmat[,i]<-x2[,i]-t(testdwdpm)
	}
	x1<-trainadjmat
	x2<-testadjmat
	cbind(x1,x2)
}

co.var <- function(x){
	100*sd(x)/mean(x)
}

SWISS <- function(data,groups,option,center){
  	d1 = nrow(data)
	n = ncol(data)
	k = max(groups)

	## Calculate group means (gmean), number of elts in each group (tempn)
  		gmeanA = matrix(nrow=d1,ncol=k)
  		tempn <-array(dim=k)
  		for(i in 1:k){
    			temp<-matrix(nrow=d1,ncol=1)
    			for(j in 1:n){
				if(groups[j]==i){
	  				temp = cbind(temp,data[,j])
				}
    			}	
  			tempn[i] = ncol(temp) - 1
  			rmeanA = rep(0,length = d1)
    			for(m in 1:d1){
      			rmeanA[m] = 0 ;
      			for (j in 1:tempn[i]){
	     				rmeanA[m] = rmeanA[m] + temp[m,j+1] ;
	  			}
      			rmeanA[m] = rmeanA[m] / tempn[i] ;
      		}
  			centdatA = matrix(nrow = d1, ncol = tempn[i])
    			for(m in 1:d1) {
      			for(j in 1:tempn[i]){
	  				centdatA[m,j] = temp[m,j+1] - rmeanA[m] 
      			}
    			}
  			gmeanA[,i] = rmeanA ;
  		}

	## Calculate overall mean (center)
		omeanA <-matrix(nrow=d1, ncol=1)
  		if(option == 1){
    			for(i in 1:d1){ 
      			omeanA[i,1] = sum(data[i,]) / n 
   			}
  		}
  		if(option == 2){
    			for(i in 1:d1){
      			omeanA[i,1] = sum(gmeanA[i,]) / k ;
    			} 
  		}
	
	## Calculate distance matrices
		distA <- matrix(nrow=2, ncol=n)
  		for(i in 1:n){
    			g = groups[i] 
    			distA[1,i] = sum((data[,i]-gmeanA[,g])^2) ;
    			distA[2,i] = sum((data[,i]-omeanA)^2) ;
   		}

	## Calculate Ratio A and B
  		RatioA = sum(distA[1,]) / sum(distA[2,]) ;
 
	#cat("SWISS") 
	#print(RatioA)
}



# SWISStest, Standardized WithIn class Sum of Squares (SWISS) Hypothesis Test
#    Permutation test for significant differences between SWISS  
#        of data1 and data2

#   Chris Cabanski's R function
# Inputs:
#     data1    - (d1 x n) matrix with d1 genes and n samples
#     data2    - (d2 x n) matrix with d2 genes and n samples
#     groups   - (n x 1) vector with labels (1,2,3,...,k) for each group
#
#     option            option for centering
#		                   1 - usual sum of squares decomposition 
#                            (uses overall mean)
#		                   2 - uses mean of the group means 
#                            (this weighs all groups equally)
#
#     nsim             Number of simulated relabellings to use
#                           (we suggest 1000)
#
# Output:
#      teststat - test statistic based on ratio of within group
#                 sum of squares divided by total sum of squares (SWICSS)
#
#      epval - empirical pvalue, based on simulated quantiles.
#                 summarizing results of permuation test
#
#      ci - 90% confidence interval of permuted population
#

#    Copyright (c) Chris Cabanski 2009

SWISStest <- function(data1,data2,groups,option,nsim){ 

  d1 = nrow(data1)
  d2 = nrow(data2)
  n = ncol(data1)
  k = max(groups)

## Calculate group means (gmean), number of elts in each group (tempn)
  gmeanA = matrix(nrow=d1,ncol=k)
  tempn <-array(dim=k)
  for(i in 1:k){
    temp<-matrix(nrow=d1,ncol=1)
    for(j in 1:n){
	if(groups[j]==i){
	  temp = cbind(temp,data1[,j])
	}
    }
  tempn[i] = ncol(temp) - 1
  rmeanA = rep(0,length = d1)
    for(m in 1:d1){
      rmeanA[m] = 0 ;
      for (j in 1:tempn[i]){
	     rmeanA[m] = rmeanA[m] + temp[m,j+1] ;
	  }
      rmeanA[m] = rmeanA[m] / tempn[i] ;
      }
  centdatA = matrix(nrow = d1, ncol = tempn[i])
    for(m in 1:d1) {
      for(j in 1:tempn[i]){
	  centdatA[m,j] = temp[m,j+1] - rmeanA[m] 
      }
    }
  gmeanA[,i] = rmeanA ;
  }

  gmeanB = matrix(nrow=d2,ncol=k)
  tempn <-array(dim=k)
  for(i in 1:k){
    temp<-matrix(nrow=d2,ncol=1)
    for(j in 1:n){
	if(groups[j]==i){
	  temp = cbind(temp,data2[,j])
	}
    }
  tempn[i] = ncol(temp) - 1
  rmeanB = rep(0,length = d2)
    for(m in 1:d2){
      rmeanB[m] = 0 ;
      for (j in 1:tempn[i]){
	     rmeanB[m] = rmeanB[m] + temp[m,j+1] ;
	  }
      rmeanB[m] = rmeanB[m] / tempn[i] ;
      }
  centdatB = matrix(nrow = d2, ncol = tempn[i])
    for(m in 1:d2){
      for(j in 1:tempn[i]){
	  centdatB[m,j] = temp[m,j+1] - rmeanB[m] 
      }
    }
  gmeanB[,i] = rmeanB ;
  }

## Calculate overall mean (center)
omeanA <-matrix(nrow=d1, ncol=1)
omeanB <-matrix(nrow=d2, ncol=1)
  if(option == 1){
    for(i in 1:d1){ 
      omeanA[i,1] = sum(data1[i,]) / n 
    }
    for(i in 1:d2){
      omeanB[i,1] = sum(data2[i,]) / n ;
    }
  }
  if(option == 2){
    for(i in 1:d1){
      omeanA[i,1] = sum(gmeanA[i,]) / k ;
    } 
    for(i in 1:d2){
      omeanB[i,1] = sum(gmeanB[i,]) / k ;
    }	
  }

## Calculate distance matrices
distA <- matrix(nrow=2, ncol=n)
distB <- matrix(nrow=2, ncol=n)
  for(i in 1:n){
    g = groups[i] 
    distA[1,i] = sum((data1[,i]-gmeanA[,g])^2) ;
    distB[1,i] = sum((data2[,i]-gmeanB[,g])^2) ;
    distA[2,i] = sum((data1[,i]-omeanA)^2) ;
    distB[2,i] = sum((data2[,i]-omeanB)^2) ;
  }


## Calculate Ratio A and B
  RatioA = sum(distA[1,]) / sum(distA[2,]) ;
  RatioB = sum(distB[1,]) / sum(distB[2,]) ;

## Standardize Matrices
   distA = distA / sum(distA[2,]) ;
   distB = distB / sum(distB[2,]) ;

## Calculate Permuted Population of Ratios

  ratio <-array(dim=nsim)

  for(i in 1:nsim ){
    perm <- runif(n)
    pdata <- matrix(nrow=2,ncol=n)
    for(j in 1:n){
      if(perm[j] < 0.5){
        pdata[,j] = distA[,j] 
      }else{ 
        pdata[,j] = distB[,j]
      }
    }
    ratio[i] = sum(pdata[1,]) / sum(pdata[2,]) ;
  }

## Calculate p-values
  countA = 0
  countB = 0
  for(i in 1:nsim){
    if(RatioA < ratio[i]){
      countA = countA + 1
    }
    if(RatioB < ratio[i]){
      countB = countB + 1
    }
  }
  epvalA = countA/nsim
  epvalA = min(epvalA,1-epvalA)
  epvalB = countB/nsim
  epvalB = min(epvalB,1-epvalB)

## Calculate 90% confidence intervals

  sorted <- sort(ratio)
  temp = floor(.05*nsim)
  left <- sorted[temp]
  temp = ceiling(.95*nsim)
  right <- sorted[temp]

  epval <- c(epvalA, epvalB)
  teststat <- c(RatioA,RatioB) 
  ci <- c(left, right) ;


  SWISS.test <- list("SWISS score for data1" = teststat[1], "SWISS score for data2" = teststat[2], "Empirical pval for data1" = epval[1], "Empirical pval for data2" = epval[2], "90% Confidence Interval" =  ci)

}

icc.2.vecs <- function (seq.vector, ma.vector)
{
	#parameters: seq.vector    - a data vector to test
	#		 ma.vector     - a data vector to compare to

	seq.n.ma  <- overlapVectors (seq.vector, ma.vector)
	seq.ratios.ma.combined<- cbind (seq.n.ma$x$dat, seq.n.ma$y$dat)
	summary(seq.ratios.ma.combined [,1])
	summary(seq.ratios.ma.combined [,2])
	icc.ratio <- icc(seq.ratios.ma.combined )$value
	rho.ratio <- cor(seq.ratios.ma.combined )[1,2]
	spearman.ratio <- cor(seq.ratios.ma.combined,method = "spearman")[1,2]
	report <- list (ICC =icc.ratio, Correlation = rho.ratio, Spearman = spearman.ratio, processed.dat = seq.ratios.ma.combined)
}


subset.by.column <- function(dm, columnHeaders)
{
#				dm		      - a data frame (nxm)
#		 		columnHeaders     - a string vector of column headers
#				return: a subset of dataframe having matched columnHeaders
#				FIXME: need to fix error when single column is requested
	headers <- dimnames(dm)[[2]]
	if (length(columnHeaders) ==1) # just single vector
	{
		index <- 0
		for (i in 1:length(headers))
		{
			if (headers[i] == columnHeaders){	index <- i	}
		}
		dt.list <- list (ID = dimnames(dm)[[1]], dat = dm[,index])
		#dt.temp <- as.data.frame(dt.list)
		header <- dimnames(dm)[[2]][index]
		#dimnames(dt.temp)[[2]] <- header
		#return (dt.temp)
		return (as.data.frame(list(dt = as.data.frame(dt.list), header = header)))
		cat ("Only single vector is selected, and usage can be different")
		str(as.data.frame(list(dt = as.data.frame(dt.list), header = header)))
	}else{
		indexes <- c()
		for (i in 1:length(columnHeaders))
		{
			indexes [i] <- header.exists(dm, columnHeaders[i])
		}
		dm.extracted <- dm[,indexes]
		return (dm.extracted)
	}
}
