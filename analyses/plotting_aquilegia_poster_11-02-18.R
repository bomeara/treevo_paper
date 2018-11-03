

setwd("d://dave//workspace//treevo_paper//")

library(ape)
library(TreEvo)



# obtain aquilegia tree (from Whittall and Hodges 2007?)
aquilegiaTree<-read.tree("datasets//aquilegia_Whttall&Hodges2007_figuredMCC.tre")
# make into a multiPhylo list
aquilegiaTreeList <- list(aquilegiaTree = aquilegiaTree)
class(aquilegiaTreeList) <- "multiPhylo"




##########################################################################################
# 2 bounds, 3 optima

# params[1] is dispersion (sigma)
# params[2] is alpha (strength of attraction to an optima)
# params[3] is rho, an exponent scaling the weighting of distance to optima
	# this parameter will control switching optima
# params[4:5] is the max boundary, for the two lower regimes regimes
# params[6:8] describes theta (optima) values for each of the three regimes

# 2 bounds, 3 optima with strong-ish attraction
params<-c(
	sigma=0.1,
	alpha=0.3,
	rho=1,
	maxbounds=c(20,40),
	theta=c(10,30,50)
	)	
	
multiOptima3IntrinsicMaxBoudary2(params=params, states=0, timefrompresent=NA)
	
		
repSim<-replicate(300,
	repeatSimSteps(params,trait = 0, nSteps = 100, 
		fun = multiOptima3IntrinsicMaxBoudary2
		)
	)
hist(repSim,main="Simulated Trait Values",breaks=20)


# same model above, with more switching between weak optima, high diffusion
params<-c(
	sigma=0.7,
	alpha=0.1,
	rho=0.5,
	maxbounds=c(20,40),
	theta=c(10,30,50)
	)	
	
multiOptima3IntrinsicMaxBoudary2(params=params, states=0, timefrompresent=NA)
	
repSim<-replicate(300,
	repeatSimSteps(params,trait = 0, nSteps = 100, 
		fun = multiOptima3IntrinsicMaxBoudary2
		)
	)
hist(repSim,main="Simulated Trait Values",breaks=20)

