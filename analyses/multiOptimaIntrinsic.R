# moser_multi-optima-single evolutionary-regime-model.R
# 08-06-18


# multi optima single evolutionary regime model
# MOSER?


#' @rdname intrinsicModels
#' @export
multiOptimaIntrinsic <- function(params, states, timefrompresent) {
    #a discrete time OU, same sd, mean, and attraction for all chars
    #params[1] is sd (sigma), params[2] is attractor (ie. character mean), params[3] is attraction (ie. alpha)
    

    sd <- params[1]
    attractor <- params[2]
    attraction <- params[3]    #in this model, this should be between zero and one
	
	#subtract current states because we want displacement
    newdisplacement <- rpgm::rpgm.rnorm(n = length(states), mean = (attractor-states)*attraction, sd = sd) 

    return(newdisplacement)
    }