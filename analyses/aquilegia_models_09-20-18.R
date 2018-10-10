
	
# September 2018
# Make two more models for Aquilegia
# 1) trait values have three optima on gradient, with some rate of switching to next-largest optima, cannot reverse
# 2) trait values evolve in three regimes with  successive upper bounds on gradient
	#(so only two upper-bounds, highest regime has no bounds) 
	# with some rate of switching to next-largest regime, cannot reverse
# addendum - maybe make rate of switching to next optima dependent on trait value?



multiOptima3IntrinsicMaxBoudary3 <- function(params, states, timefrompresent) {
	# 1) trait values have three optima on gradient, with some rate of
		# switching to next-largest optima, cannot reverse
	#
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
	# BUT this also has max bounds!
	# proximity to theta determines which theta a lineage is in, but lineages 
	#
	sigma<-params[1]
	alpha<-params[2]
	rho <- params[3]
	maxBound<-params[4:6][i]
	theta<-params[7:9][i]	
	#
	# what regime does the lineage sit in?
	# (1) cannot be in a regime whose max bound it has surpassed
	unsurpassedBound <- maxBound >= states	
	# also make it so that lineages in the first regime can't just to the third
	if(all(unsurpassedBound)){
		unsurpassedBound[3] <- FALSE
		}
	#
	# (2) it has a chance of being in the next highest regime 
	# calculate probabilistic weights relative to distance from all theta
		# raised to the power of rho - scaling parameter
	thetaWeights <- (1/abs(theta-states))^rho
	#
	# (3) combined 1+2, rescale so sum to 1 as probability
	thetaWeights <- thetaWeights*as.numeric(unsurpassedBound)
	thetaWeights<-thetaWeights/sum(thetaWeights)
	#
	# sample a theta/bound regime
	regime<-sample(1:length(theta), 1, prob = thetaWeights)
	theta<-theta[regime]
	maxBound<-maxBound[regime]
	#
	# now 
	#subtract current states because we want displacement
    newdisplacement <- rpgm::rpgm.rnorm(
		n = length(states),
		mean = (theta-states)*alpha,
		sd = sd) 
	# shift states if they surpass max bound
	for (i in 1:length(newdisplacement)) {
        newstate <- newdisplacement[i]+states[i]
        if (newstate>maxBound) { #newstate less than min
            newdisplacement[i] <- maxBound-states[i] 
			#so, rather than go above the maximum, this moves the new state to the maximum
			}
		}
	#
    return(newdisplacement)
    }


multiOptima3IntrinsicMaxBoudary2 <- function(params, states, timefrompresent) {
	# trait values evolve in three regimes with successive upper bounds on gradient (
	   # so only ****two**** upper-bounds, highest regime has no bounds) 
		# with some rate of switching to next-largest regime, cannot reverse
	#
    #a discrete time OU with multiple optima in three regimes
		# with equal attraction (alpha) to all optima (theta 1:N)
		# and each regime having its own max trait value
	# breakdown of params:
		# params[1] is dispersion (sigma)
		# params[2] is alpha (strength of attraction to an optima)
		# params[3] is rho, an exponent scaling the weighting of distance to optima
			# this parameter will control switching optima
		# params[4:5] is the max boundary, for the two lower regimes regimes
		# params[7:9] describes theta (optima) values for each of the three regimes
	# In this model, optima represent fixed trait values conveying adaptive benefit
		# the proximity of a population to an optima makes it more likely to be under that regime
		# a point equidistant between multiple regimes may be drawn to any
	# the draw to any specific optima is inverse to distance from optima
	# thus a lineage at an optima may show large variance as it circles the plateau
		# then suddenly feel drawn to another optima, and show sudden, giant shifts toward that optima
	# BUT this also has max bounds!
	# proximity to theta determines which theta a lineage is in, but lineages 
	#
	sigma<-params[1]
	alpha<-params[2]
	rho <- params[3]
	# add an infinite max bound to the last regime
	maxBound<-c(params[4:5],Inf)
	theta<-params[7:9]	
	#
	# what regime does the lineage sit in?
	# (1) cannot be in a regime whose max bound it has surpassed
	unsurpassedBound <- maxBound >= states	
	# also make it so that lineages in the first regime can't just to the third
	if(all(unsurpassedBound)){
		unsurpassedBound[3] <- FALSE
		}
	#
	# (2) it has a chance of being in the next highest regime 
	# calculate probabilistic weights relative to distance from all theta
		# raised to the power of rho - scaling parameter
	thetaWeights <- (1/abs(theta-states))^rho
	#
	# (3) combined 1+2, rescale so sum to 1 as probability
	thetaWeights <- thetaWeights*as.numeric(unsurpassedBound)
	thetaWeights<-thetaWeights/sum(thetaWeights)
	#
	# sample a theta/bound regime
	regime<-sample(1:length(theta), 1, prob = thetaWeights)
	theta<-theta[regime]
	#
	# now 
	#subtract current states because we want displacement
    newdisplacement <- rpgm::rpgm.rnorm(
		n = length(states),
		mean = (theta-states)*alpha,
		sd = sd) 
	# shift states if they surpass max bound
	if(regime != 3){
		maxBound<-maxBound[regime]
		for (i in 1:length(newdisplacement)) {
			newstate <- newdisplacement[i]+states[i]
			if (newstate>maxBound) { #newstate less than min
				newdisplacement[i] <- maxBound-states[i] 
				#so, rather than go above the maximum, this moves the new state to the maximum
				}
			}
		}
	#
    return(newdisplacement)
    }


	
	
	
# old aquilegia model from O'Meara, slightly rewritten
	# levy walk model
pollinatorShiftIntrinsic <-function(params, states, timefrompresent) {
  #params[1] is sd for normal process, 
  #params[2] is sd for jump size of non-background shift 
  #params[3] is mean of shift (mean of normal process is zero), 
  #params[4] is per generation probability of a non-background shift/jump
  if (runif(1,0,1) < params[4] ) {
    # *not* background.process
	newdisplacement<-rnorm(n=length(states),mean=params[3],sd=params[2])
    return(newdisplacement)
  }else{
	# background process
    newdisplacement<-rnorm(n=length(states),mean=0,sd=params[1])
    return(newdisplacement)
	}
  }
  

# designating priors for old aquilegia model from legacy code
# from B O'Meara
# 
#do fixed for param2 b/c dont have a lot of data to use
# intrinsicPriorsFns=c("exponential","fixed","uniform","uniform") 
#
#intrinsicPriorsValues=matrix(
#	c(
#		# regular rate/sd of change = 0.08 is based on best guess of average rate of change, it seems
#		rep(1/0.08 , 2),
#		# jump rate of change - fixed at zero because:
#			#"do fixed for param2 b/c dont have a lot of data to use" ??
#		c(0,0),
#		# mean of shift (if not zero, then its a trend model)
#		c(-maxCharDistance, maxCharDistance),
#		# rate of a shift
#		0, 0.00001
#		)
#	, nrow=2, byrow=FALSE)
#
# hmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmm no

