rm(list=ls())
library(TreEvo)
library(ape)
library(geiger)
library(parallel)
ntax<-30
abcMethod<-"rejection"
abcTolerance<-0.05

birthdeath.tree(1,0, taxa.stop=ntax+1)->l
drop.tip(l, as.character(ntax+1))->phy
phy$edge.length<-phy$edge.length/max(branching.times(phy))

intrinsicFn=brownianIntrinsic
extrinsicFn=nullExtrinsic
startingPriorsFns="uniform"
startingPriorsValues=matrix(c(0,60),nrow=2)
intrinsicPriorsFns=c("uniform")
intrinsicPriorsValues=matrix(c(0, .10), nrow=2, byrow=FALSE)
extrinsicPriorsFns=c("fixed")
extrinsicPriorsValues=matrix(c(0, 0), nrow=2, byrow=FALSE)

startingValuesTrue<-c(10)
intrinsicValuesTrue<-c(0.05)
extrinsicValuesTrue<-c(0)
TreeYears=1e+04
timeStep<-1/TreeYears
#nloops<-10000
#nsimsPerLoop<-12*6
nrepSims<-10000000
all.results<-list()



splits<-getSimulationSplits(phy)

char<-convertTaxonFrameToGeigerData(
	doSimulation(
		splits=splits, 
		intrinsicFn=intrinsicFn,
		extrinsicFn=extrinsicFn,
		startingValues=startingValuesTrue, 
		intrinsicValues=intrinsicValuesTrue,
		extrinsicValues=extrinsicValuesTrue,
		timeStep=timeStep,
		saveRealParams=T,
		jobName="push"), phy
)


#for (loopID in sequence(nloops)) {
	trial2<-doRun_rej(
		phy=phy, 
		traits=char,
		intrinsicFn=intrinsicFn, 
		extrinsicFn=extrinsicFn, 
		startingPriorsFns=startingPriorsFns,
		startingPriorsValues=startingPriorsValues,
		intrinsicPriorsFns=intrinsicPriorsFns,
		intrinsicPriorsValues=intrinsicPriorsValues,
		extrinsicPriorsFns=extrinsicPriorsFns,
		extrinsicPriorsValues=extrinsicPriorsValues,
		StartSims=nrepSims, 
		jobName="trial2", 
		abcTolerance=abcTolerance, 
		multicore=T, 
		coreLimit=detectCores(),
		TreeYears=1/timeStep,
		checkpointFile="ntax30_bm_rej",
		checkpointFreq=5*detectCores()
	)	
#	all.results[[length(all.results)+1]]<-trial2
#	print(paste("length all.results is",length(all.results)))
#	print(paste("number of completed sims is",nsimsPerLoop*loopID))
#	warnings()
#	try(system("mv ntax30_bm_rej.RSave ntax30_bm_rej.previous.RSave"))
	save(list=ls(),file="ntax30_bm_rej.rda")
#}	


