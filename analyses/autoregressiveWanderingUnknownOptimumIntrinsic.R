#' @rdname intrinsicModels
#' @export
autoregressiveWanderingUnknownOptimumIntrinsic <- function(params, states, timefrompresent) {
	# 3) Time - autoregressive model with optimum based on a
		# factor that changes through time (like O2 concentration)
	# Parameters of the regression: 
		# a single variable function to convert O2 to optimal gill size or whatever
		# strength of pull
		# BM wiggle
	
	# 10-20-18
	# can't do it with an unknown optimum because the model function doesn't talk to other
	# instances of the model, so how would they know what the optimum is at a particular time-point
	# ->>>>>>>>>>>>>>>> 
	#       We need a dataset with an environmental dataset
	#            that we can treat as an optimum being tracked
	# but what? and this will require another argument.
		



    #a discrete time OU, same sd, mean, and attraction for all chars
    #params[1] is sd (sigma), params[2] is attractor (ie. character mean), params[3] is attraction (ie. alpha)
    sd <- params[1]
    attractor <- params[2]
    attraction <- params[3]    #in this model, this should be between zero and one
    #subtract current states because we want displacement
	newdisplacement <- rpgm::rpgm.rnorm(
					n = length(states), 
					mean = (attractor-states)*attraction, 
					sd = sd) 
    return(newdisplacement)
} 