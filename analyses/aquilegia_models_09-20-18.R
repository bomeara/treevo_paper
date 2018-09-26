

multiOptimaIntrinsicMaxBoudary3 <- function(params, states, timefrompresent) {
    #a discrete time OU with multiple optima in three regimes
		# with equal attraction (alpha) to all optima (theta 1:N)
		# and each regime having its own max trait value
	# breakdown of params:
		# params[1] is dispersion (sigma)
		# params[2] is alpha (strength of attraction to an optima)
		# params[3] is rho, an exponent scaling the weighting of distance to optima
			# this parameter will control switching optima
		# params[4:6] is the max boundary, for each of the three regimes
		# params[7:9] describes theta (optima) values for each of the three regimes
		
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
	
maxBoundaryAutoregressiveIntrinsic <- function(params, states, timefrompresent) {
    #a discrete time OU, same sd, mean, and attraction for all chars
    #params[1] is sd (sigma), 
		# params[2] is attractor (ie. character mean), params[3] is attraction (ie. alpha), 
		# params[4] is max bound
    sd <- params[1]
    attractor <- params[2]
    attraction <- params[3]    #in this model, this should be between zero and one
    minBound <- params[4]
    newdisplacement <- rpgm::rpgm.rnorm(
		n = length(states), 
		mean = (attractor-states)*attraction, 
		sd = sd) #subtract current states because we want displacement
    #message(newdisplacement)
	for (i in length(newdisplacement)) {
        newstate <- newdisplacement[i]+states[i]
        if (newstate>params[2]) { #newstate less than min
            newdisplacement[i] <- params[2]-states[i] 
			#so, rather than go below the minimum, this moves the new state to the maximum
        }
    }
    return(newdisplacement)
}



	
    #a discrete time OU, same sd, mean, and attraction for all chars
    #params[1] is sd (sigma), 
		# params[2] is attractor (ie. character mean), params[3] is attraction (ie. alpha), 
		# params[4] is max bound	
	#
	# September 2018
	# Make two more models for Aquilegia
	# 1) trait values have three optima on gradient, with some rate of switching to next-largest optima, cannot reverse
	# 2) trait values evolve in three regimes with  successive upper bounds on gradient (so only two upper-bounds, highest regime has no bounds) with some rate of switching to next-largest regime, cannot reverse
	# addendum - maybe make rate of switching to next optima dependent on trait value?

	
	
	
	# September 2018
	# Make two more models for Aquilegia
	# 1) trait values have three optima on gradient, with some rate of switching to next-largest optima, cannot reverse
	# 2) trait values evolve in three regimes with  successive upper bounds on gradient (so only two upper-bounds, highest regime has no bounds) with some rate of switching to next-largest regime, cannot reverse
	# addendum - maybe make rate of switching to next optima dependent on trait value?