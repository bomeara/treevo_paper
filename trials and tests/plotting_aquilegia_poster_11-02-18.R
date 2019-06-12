
# aquilegia_models_test_10-10-18.R

library(ape)
library(TreEvo)

source("D:\\dave\\workspace\\treevo_paper\\analyses\\aquilegia_models_09-20-18.R")

# simulate n time-steps, repeat many times, plot results
repeatSimSteps<-function(params,trait=0,nSteps,fun){
	for(i in 1:nSteps){
	# add to original trait value to get new trait value
		trait<-trait+fun(
			params=params, states=trait, timefrompresent=NA)
		}
	trait
	}

set.seed(1)

setwd("d://dave//workspace//treevo_paper//")


# obtain aquilegia tree (from Whittall and Hodges 2007?)
aquilegiaTree<-read.tree("datasets//aquilegia_Whttall&Hodges2007_figuredMCC.tre")
# make into a multiPhylo list
aquilegiaTreeList <- list(aquilegiaTree = aquilegiaTree)
class(aquilegiaTreeList) <- "multiPhylo"


# obtain aquilegia trait data (from Whittall and Hodges 2007?)
	# need both nectur spur lengths and regime data
	#
aquilegiaTrait<-read.table("datasets//aquilegia_traitData.txt", header=FALSE, row.names=1)

# get just nectur spur length
aquilegiaSpurLength<-aquilegiaTrait[,2]
# and take the natural log
	# (note that the third column of the table was already the natural log)
	# previous code from Brian had 'log(data[,3])' - log of a log
aquilegiaSpurLength<-log(aquilegiaSpurLength)


# aquilegia regimes - pollinator syndromes
aquilegiaPollinators<-aquilegiaTrait[,14]
# regimes coded 0, 1, 2
	# 0 is bumble-bee, 1 is humming-bird, 2 is hawkmoth
# this probably won't be used directly?
# could use for post-analysis comparisons? Hmm

library(ggplot2)
aqData<-data.frame(aquilegiaSpurLength=aquilegiaSpurLength,
	aquilegiaPollinators=as.factor(aquilegiaPollinators)
	)
ggplot(aqData,
   aes(x=aquilegiaSpurLength, fill=aquilegiaPollinators)) +
  geom_histogram(show.legend=FALSE)


#plot the entire data set (everything)
breaks<-seq(1,5,by=0.5)

hist(aquilegiaSpurLength, breaks=breaks, col="Yellow",main="")

#then everything except one sub group (1 in this case)
hist(aquilegiaSpurLength[aquilegiaPollinators!=2],
	 breaks=breaks, col="Red", add=TRUE)

#then everything except two sub groups (1&2 in this case)
hist(aquilegiaSpurLength[aquilegiaPollinators!=2 & aquilegiaPollinators!=1],
	 breaks=breaks, col="Blue", add=TRUE)



##########################################################################################
# 2 bounds, 3 optima

# params[1] is dispersion (sigma)
# params[2] is alpha (strength of attraction to an optima)
# params[3] is rho, an exponent scaling the weighting of distance to optima
	# this parameter will control switching optima
# params[4:5] is the max boundary, for the two lower regimes regimes
# params[6:8] describes theta (optima) values for each of the three regimes


layout(matrix(1:12,4,3))
par(mar=c(2,2,0,0))

for(i in 1:12){
# same model above, with more switching between weak optima, high diffusion
params<-c(
	sigma=runif(1,0.2,0.4),
	alpha=0.1,
	rho=runif(1,0.6,0.8),
	maxbounds=c(20,40),
	theta=c(10,30,50)
	)	
	
repSim<-replicate(31,
	repeatSimSteps(params,trait = 0, nSteps = 100, 
		fun = multiOptima3IntrinsicMaxBoundary2
		)
	)
hist(repSim,main="Simulated Trait Values",
	breaks=20,axes=FALSE,ann=FALSE)
Axis(side=1, labels=FALSE)
Axis(side=2, labels=FALSE)

}


