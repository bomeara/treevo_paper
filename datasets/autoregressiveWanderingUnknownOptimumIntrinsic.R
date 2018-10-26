#' @rdname intrinsicModels
#' @export
autoregressiveWanderingUnknownOptimumIntrinsic <- function(params, states, timefrompresent) {
	# 3) Time - autoregressive model with optimum based on a
		# factor that changes through time (like O2 concentration)
	# and presumably the tracked env factor will be analyzed many times over
	#
	# Parameters of the regression: 
		# a single variable function to convert O2 to optimal gill size or whatever
		# strength of pull
		# BM wiggle
	#
	# 10-20-18
	# can't do it with an unknown optimum because the model function doesn't talk to other
	# instances of the model, so how would they know what the optimum is at a particular time-point
	# ->>>>>>>>>>>>>>>> 
	#       We need a dataset with an environmental dataset
	#            that we can treat as an optimum being tracked
	# but what data? this will require another argument.
	#	
	# 10-23-18
	# Previous discussions of what models we wanted to test included discussion a "time / autoregressive" model
	# as far as I can tell from our previous conversations, we had meant for this to be an autoregressive model where the optimum follows some known environmental predictor through time
	# (note that the way treevo is designed, it is impossible as far as I can tell to do an arbitrary environmental predictor)
	# (i.e. an unknown autoregressive optima that varies over time)
	# a model that follows a known environmental predictor could be coded
	# it would require some small but not minor changes to the treevo code base, particularly that the variable would need to be an extra argument to all intrinsic models
	# and the way time is passed to the intrinsic models would need to be in user-defined time units - i.e. the same units as the tree / matrix of enviromental variable itself
	# Previous discussions of what models we wanted to test included discussion a "time / autoregressive" model
	# as far as I can tell from our previous conversations, we had meant for this to be an autoregressive model where the optimum follows some known environmental predictor through time
	# (note that the way treevo is designed, it is impossible as far as I can tell to do an arbitrary environmental predictor)
	# (i.e. an unknown autoregressive optima that varies over time)
	# a model that follows a known environmental predictor could be coded, but it would require some small but not minor changes to the treevo code base, particularly that the variable would need to be an extra argument to all intrinsic models
	# and the way time is passed to the intrinsic models would need to be in user-defined time units - i.e. the same units as the tree / matrix of enviromental variable itself
	#
	#
	# Brian 10-23-18
	#
	# environVariable DEFINED within the function - then doesn't need to be input carried through doRun
		# and instead is in the function by default
		# okay, that should work
	
	
	
	
	# need to get the environmental variable
	environVariable[timefrompresent]


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
