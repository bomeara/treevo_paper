# aquilegia_models_test_10-10-18.R

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

#########################################################################
# 3 bounds, 3 optima

# params[1] is dispersion (sigma)
# params[2] is alpha (strength of attraction to an optima)
# params[3] is rho, an exponent scaling the weighting of distance to optima
	# this parameter will control switching optima
# params[4:6] is the max boundary, for each of the three regimes
# params[7:9] describes theta (optima) values for each of the three regimes

# 3 bounds, 3 optima with weak-ish attraction
params<-c(
	sigma=0.1,
	alpha=0.3,
	rho=1,
	maxbounds=c(20,40,60),
	theta=c(10,30,50)
	)	
	
multiOptima3IntrinsicMaxBoudary3(params=params, states=0, timefrompresent=NA)
	
	
repSim<-replicate(300,
	repeatSimSteps(params,trait = 0, nSteps = 100, 
		fun = multiOptima3IntrinsicMaxBoudary3
		)
	)
hist(repSim,main="Simulated Trait Values",breaks=20)


# same model above, with more switching between optima
params<-c(
	sigma=0.1,
	alpha=0.5,
	rho=0.3,
	maxbounds=c(20,40,60),
	theta=c(10,30,50)
	)	
	
multiOptima3IntrinsicMaxBoudary3(params=params, states=0, timefrompresent=NA)
	
repSim<-replicate(300,
	repeatSimSteps(params,trait = 0, nSteps = 100, 
		fun = multiOptima3IntrinsicMaxBoudary3
		)
	)
hist(repSim,main="Simulated Trait Values",breaks=20)


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

