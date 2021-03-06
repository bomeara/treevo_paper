# moser_multi-optima-single evolutionary-regime-model.R
# 08-06-18


# multi optima single evolutionary regime model
# MOSER?


#' @rdname intrinsicModels
#' @export
multiOptimaIntrinsic <- function(params, states, timefrompresent) {
    #a discrete time OU with multiple optima in the same regime 
		# with equal attraction (alpha) to all optima (theta 1:N)
	# breakdown of params:
		# params[1] is dispersion (sigma)
		# params[2] is alpha (strength of attraction to an optima)
		# params[3] is rho, an exponent scaling the weighting of distance to optima
			# this parameter will control switching optima
		# params[4:n] describes theta values
			# n-2 = N # of optima describe by this model
	# In this model, optima represent fixed trait values conveying adaptive benefit
		# the proximity of a population to an optima makes it more likely to be under that regime
		# a point equidistant between multiple regimes may be drawn to any
	# the draw to any specific optima is inverse to distance from optima
	# thus a lineage at an optima may show large variance as it circles the plateau
		# then suddenly feel drawn to another optima, and show sudden, giant shifts toward that optima
	# this all seems realistic...
	#
	sigma<-params[1]
	alpha<-params[2]
	rho <- params[3]
	theta<-params[-(1:3)]
	#
	# measure distances to theta
	# convert to probabilistic weights
		# raised to the power of rho - scaling parameter
	thetaWeights<-(1/abs(theta-states))^rho
	# rescale so sum to 1, as probabilities
	thetaWeights<-thetaWeights/sum(thetaWeights)
	# sample a theta
	theta<-sample(theta,1,prob=thetaWeights)
	# now 
	#subtract current states because we want displacement
    newdisplacement <- rpgm::rpgm.rnorm(n = length(states), mean = (theta-states)*alpha, sd = sd) 
    return(newdisplacement)
    }

# three optima model, with strong attraction	
set.seed(1)
params<-c(
	sigma=0.1,
	alpha=0.7,
	rho=1,
	theta=c(-20,20,50)
	)	
	
multiOptimaIntrinsic(params=params, states=0, timefrompresent=NA)

# simulate n time-steps, repeat many times, plot results
repeatSimSteps<-function(params,trait=0,nSteps){
	for(i in 1:nSteps){
	# add to original trait value to get new trait value
		trait<-trait+multiOptimaIntrinsic(
			params=params, states=trait, timefrompresent=NA)
		}
	trait
	}
repSim<-replicate(300,repeatSimSteps(params,trait=0,100))
hist(repSim,main="Simulated Trait Values",breaks=20)


# same model above, with more switching between optima
set.seed(1)
params<-c(
	sigma=0.1,
	alpha=0.7,
	rho=0.5,
	theta=c(-20,20,50)
	)	
	
multiOptimaIntrinsic(params=params, states=0, timefrompresent=NA)

# simulate n time-steps, repeat many times, plot results
repeatSimSteps<-function(params,trait=0,nSteps){
	for(i in 1:nSteps){
	# add to original trait value to get new trait value
		trait<-trait+multiOptimaIntrinsic(
			params=params, states=trait, timefrompresent=NA)
		}
	trait
	}
repSim<-replicate(300,repeatSimSteps(params,trait=0,100))
hist(repSim,main="Simulated Trait Values",breaks=20)
