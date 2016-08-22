rm(list=ls())
try(detach("package:TreEvo",unload=TRUE))
install.packages("/Users/bomeara/Documents/MyDocuments/Active/TreEvo/pkg",repos=NULL,type="source")
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



SVtoPlot<-colnames(all.SV)[c(1:23,53,82,111)]

for (tFV.index in sequence(dim(all.tFV)[2])) {
  pdf(file=paste("raw.SV",tFV.index,"pdf",sep="."),width=15,height=10)
  par(mfrow=c(4,7))
  for (SV.index in sequence(length(SVtoPlot))) {
    SV.col<-which(colnames(all.boxcox.SV) == SVtoPlot[SV.index])
    col.main="darkgray"
    if(whichVip[SV.col,tFV.index]) {
      col.main="red" 
    }
    plot(x=all.tFV[,tFV.index],y=abs(all.boxcox.SV[,SV.col]-original.boxcox.SV[SV.col]),xlab="",ylab="",bty="n",main=paste(SVtoPlot[SV.index]," (",round(allVip[SV.col,tFV.index],2),")",sep=""),pch='.',col=rgb(0,0,0,0.1),col.main=col.main)
    lines(rep(original.tFV[tFV.index],2),range(abs(all.boxcox.SV[,SV.col]-original.boxcox.SV[SV.col])),col="red")
  }
  dev.off()
  system(paste("open ",paste("raw.SV",tFV.index,"pdf",sep=".")))
}

#print(cbind(original.SV,apply(all.SV,2,mean),original.boxcox.SV,apply(all.boxcox.SV,2,mean)))


trueFreeValuesMatrix<-all.tFV
boxcoxSummaryValuesMatrix<-all.boxcox.SV
boxcoxOriginalSummaryStats<-original.boxcox.SV

abcDistancesRaw<-matrix(nrow=dim(trueFreeValuesMatrix)[1], ncol=dim(trueFreeValuesMatrix)[2])  #rep(0,dim(trueFreeValuesMatrix)[1]) #will hold the distances for each particle
#abcDistancesRawTotal<-vector(length=dim(trueFreeValuesMatrix)[1])
#now we go true parameter by true parameter, using only the summary values with enough importance for each
#we get a distance for each particle from the observed particle for each true param
#then get a euclidean distance for all of these
#it's like getting delta X, then delta Y, and figuring out distance from the origin using
#sqrt((X-0)^2 + (Y-0)^2)
for (freeParamIndex in sequence(dim(trueFreeValuesMatrix)[2])) {
  abcDistancesRaw[,freeParamIndex]<-abc(target=boxcoxOriginalSummaryStats[whichVip[,freeParamIndex]], param=trueFreeValuesMatrix[,freeParamIndex], sumstat= boxcoxSummaryValuesMatrix[,whichVip[,freeParamIndex]], tol=1, method=abcMethod)$dist^2 #because we're doing euclidean distance, from observed, which has value 0, 0, 0, etc.
}
abcDistancesRawTotal<-apply(abcDistancesRaw, 1, sum)
abcDistances<-sqrt(abcDistancesRawTotal) #Euclid rules.

abcResults<-vector("list")
abcResults$unadj.values<-trueFreeValuesMatrix[which(abcDistances<=quantile(abcDistances, prob=abcTolerance)), ] #here's where we diy abc
abcResults$dist<-abcDistances[which(abcDistances<=quantile(abcDistances, prob=abcTolerance))]

#special<-round(runif(1,min=1,max=dim(all.SV)[1]))
#focal.value<-all.SV[special,1]
#focal.value<-original.SV[1]
#true.value<-all.tFV[special,2]
#true.value<-original.tFV[2]
#plot(all.tFV[,2],abs(focal.value-all.SV[,1]))
#lines(rep(true.value,2),range(abs(focal.value-all.SV[,1])),col="red")

#original.SV[1]
#range(all.SV[,1])
#plot(density(all.SV[,1]))
#plot(all.SV[,1],all.SV[,2])



pdf(file="RootState.pdf",width=10,height=5)
priorCurve=getUnivariatePriorCurve(priorValues=c(0,60), priorFn="uniform")
posteriorCurve=getUnivariatePosteriorCurve(abcResults$unadj.values[,1],from=min(priorCurve$x),to=max(priorCurve$x))
plotUnivariatePosteriorVsPrior(posteriorCurve, priorCurve, label="Root state", trueValue=10, prob=0.95)
dev.off()
system("open RootState.pdf")

pdf(file="BMrate.pdf",width=10,height=5)
priorCurve=getUnivariatePriorCurve(priorValues=c(0,0.1), priorFn="uniform")
posteriorCurve=getUnivariatePosteriorCurve(abcResults$unadj.values[,2],from=min(priorCurve$x),to=max(priorCurve$x))
plotUnivariatePosteriorVsPrior(posteriorCurve, priorCurve, label="Brownian rate", trueValue=0.05, prob=0.95)
dev.off()
system("open BMrate.pdf")
