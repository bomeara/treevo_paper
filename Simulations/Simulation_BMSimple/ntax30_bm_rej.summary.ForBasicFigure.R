rm(list=ls())
#try(detach("package:TreEvo",unload=TRUE))
#install.packages("/Users/bomeara/Documents/MyDocuments/Active/TreEvo/pkg",repos=NULL,type="source")
source("/Users/bomeara/Documents/MyDocuments/Active/TreEvo/pkg/R/plotUnivariatePosteriorVsPrior.R")

library(TreEvo)
setwd("/Users/bomeara/Documents/MyDocuments/Ongoing/Talks/2012vii07_Evolution2012_TreEvo/Simulation_BMSimple2_runL")
load("ntax30_bm_rej.intialsettings.Rsave")
load("ntax30_bm_rej.trueFreeValuesANDSummaryValues.Rsave")

#Note that these are copied from the batching function; verify that they match
startingValuesTrue<-c(10)
intrinsicValuesTrue<-c(0.05)
extrinsicValuesTrue<-c(0)
TreeYears=1e+04
timeStep<-1/TreeYears

vipthresh=0.8
abcMethod="rejection"
abcTolerance=0.01
print(dim(trueFreeValuesANDSummaryValues))
all.tFV<-trueFreeValuesANDSummaryValues[,1:2]
all.SV<-trueFreeValuesANDSummaryValues[,c(-1,-2)]
colnames(all.SV)<-sumStatNames(phy)

boxcoxEstimates<-boxcoxEstimation(all.SV)
boxcoxAddition<-boxcoxEstimates$boxcoxAddition
boxcoxLambda<-boxcoxEstimates$boxcoxLambda
all.boxcox.SV<-boxcoxEstimates$boxcoxSummaryValuesMatrix
colnames(all.boxcox.SV)<-sumStatNames(phy)

char<-traits
original.SV<-summaryStatsLong(phy=phy, data=char)
names(original.SV)<-sumStatNames(phy)
original.boxcox.SV<-boxcoxTransformation(original.SV, boxcoxAddition, boxcoxLambda)
names(original.boxcox.SV)<-sumStatNames(phy)
original.tFV<-c(startingValuesTrue,intrinsicValuesTrue)

library("mixOmics")
getVipSingleColumn<-function(trueFreeValuesColumn, boxcoxSummaryValuesMatrix) {
  return(vip(pls(X=boxcoxSummaryValuesMatrix, Y=trueFreeValuesColumn, ncomp=1)))
}

#note this has vip for summary val 5, for true param 1, in row 5, column 1, and so forth
allVip<-apply(all.tFV, 2, getVipSingleColumn, all.boxcox.SV)

#this will have which summary vals have importance above the threshold
whichVip<-(allVip>vipthresh)


plotQuantile<-function(p=0.5) {
special<-which.min(abs(all.tFV[,2]-quantile(all.tFV[,2],p)))
original.SV<-all.SV[special,]
#names(original.SV)<-sumStatNames(phy)
original.boxcox.SV<-boxcoxTransformation(original.SV, boxcoxAddition, boxcoxLambda)
#names(original.boxcox.SV)<-sumStatNames(phy)
original.tFV<-all.tFV[special,]


#SVtoPlot<-colnames(all.SV)[c(1:23,53,82,111)]
SVtoPlot<-colnames(all.SV)[11]
tFV.index<-2
  #pdf(file=paste("basic.SV",tFV.index,"pdf",sep="."),width=15,height=10)
  #  par(mfrow=c(4,7))
  for (SV.index in sequence(length(SVtoPlot))) {
    SV.col<-4
    col.main="darkgray"
    if(whichVip[SV.col,tFV.index]) {
      col.main="red" 
    }
    plot(x=all.tFV[,tFV.index],y=abs(all.boxcox.SV[,SV.col]-original.boxcox.SV[SV.col]),xlab="",ylab="",bty="n",pch='.',col=rgb(0,0,0,0.4),yaxt="n")
    lines(rep(original.tFV[tFV.index],2),range(abs(all.boxcox.SV[,SV.col]-original.boxcox.SV[SV.col])),col="red")
  }
  #dev.off()
  #system(paste("open ",paste("basic.SV",tFV.index,"pdf",sep=".")))
}

#print(cbind(original.SV,apply(all.SV,2,mean),original.boxcox.SV,apply(all.boxcox.SV,2,mean)))
