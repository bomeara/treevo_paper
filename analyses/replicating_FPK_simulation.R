

# actualistic bounds on function
getTraitBounds<-function(x){
	bounds <- c(
		min(x)-((max(x) - min(x))/2),
		max(x)+((max(x) - min(x))/2)
		)
	return(bounds)
	}


traitData<-rnorm(100,0,1)
# need traits to calculate bounds
bounds<-getTraitBounds(traitData)

# example from vignette for 
	# should make two peak landscape
params<-c(
	a=1,
	b=0.5,
	c=0,
	sigma=1,
	bounds)
trait<-traitData[1]


#' @rdname intrinsicModels
#' @export
landscapeFPK_Intrinsic <- function(params, states, timefrompresent) {
	#a discrete time FPK


# From Boucher et al:
# Finally, note that both BM and the OU model are special cases of the FPK
# model: BM corresponds to V(x)=0 and OU to
# V(x)=((alpha/sigma^2)*x^2)-((2*alpha*theta/(sigma^2))*x)

	########################################################
	# equation for getting potential	
	potentialFun<-function(x,a,b,c){
		# V(x)=ax4+bx2+cx 
		Vres <- (a*(x^4))+(b*(x^2))+(c*x)
		return(Vres)
		}

	############################################

	# parameters: a,b,c, rootState, sigma
		# unnecc: rootState 
	sigma <- params[4]
	bounds <- params[5:6]

	# sim controls
	grainscaleFPK <- 10
	
grainScale<-grainscaleFPK
# rename M
# rename eigM



##############################################################################

# all of the following only need to be run
	# when the parameters of FPK are changed
# this can be pre-calculated for a single run with lexical scoping	
	
# landscape descriptor function
# over the arbitrary interval (-1.5 : 1.5)
arbSequence<-seq(from=-1.5,to=1.5,
	length.out=grainScale)
origSequence<-seq(from=bounds[1],to=bounds[2],
	length.out=grainScale)
# # potentialVector is numeric vector representing the potential
	# length = grainScale
# V(x)=ax4+bx2+cx 
potentialVector<-potentialFun(
	x=arbSequence,
	a=params[1],b=params[2],c=params[3])	
#
# Coefficient of Diffusion of the Model
dCoeff <- log((sigma)^2/2) 
#
# Transition matrix describing prob of evolving between two sites in the trait grid in an infinitesimal time step.	
# Create and diagonalize the transition matrix that has been discretized
# returns: the transition matrix going forward in time, for simulating traits only
#
# make empty matrix
expD <- M <- matrix(0,grainScale,grainScale)
#assign values not on outer rows/columns
for (i in 1:(grainScale)){
	if(i>1){
		M[i-1,i] <- exp((potentialVector[i]-potentialVector[i-1])/2)
		}
	if(i<grainScale){
		M[i+1,i] <- exp((potentialVector[i]-potentialVector[i+1])/2)
		}
	# rate of staying in place is negative sum of neighboring cells
	neighbors<-c(ifelse(i>1,M[i-1,i],0),
			ifelse(i<grainScale,M[i+1,i],0)
			)
	M[i,i] <- (-sum(neighbors))
	}

# eigenvalues and eigenvectors of transition matrix
	# take only real components
eigM <- lapply(eigen(M),Re)
# 			
# solve the eigenvectors
solvedEigenvectors <- solve(eigM$vectors,tol = 1e-30) 
#
# get expected dispersion
# scale expected dispersion to original trait scale
	# (tau from Boucher et al.'s original code)
origScaler <- (((bounds[2]-bounds[1])/(grainScale-1))^2)
# assign dispersion to diagonal of expD
diag(expD) <- exp(exp(dCoeff)/origScaler*eigM$values)
# previous time-dep version from Boucher et al's code
	# diag(expD) <- exp(t*diag_expD)
#
# take dot product of expD, eigenvectors and solved eigenvectors
# get matrix of potential for future trait values
	# given potentialVector alone, not yet considering previous values
potentialMatrix <- eigM$vectors%*%expD%*%solvedEigenvectors

###############################################



###############################################

	x <- trait
	

#t <- tree$edge.length[i]


	origIntLength<-(grainScale-1)/(bounds[2]-bounds[1])
	

	# need a vector, length = grainscale
		# with zeroes in intervals far from current trait value
		# and with density=1 distributed in interval=origIntLength
			# around the original trait value


	intDensity <- rep(0,grainScale)  
	traitInterval<-c(trait-origIntLength/2, trait+origIntLength/2)
	whichInts<-sapply(traitInterval,function(x) min(which(origSequence>=x),grainScale))

	#if(length(whichInts)!=2){
	#	stop("something that is length of one interval should be in two intervals")
	#	}

	origSequence[whichInts[1]]+traitInterval
	intDensity<-ifelse(trait>origSequence[-1],
		origSequence[-1]-traitInterval[1],
		traitInterval[2]-origSequence[-grainScale]
		)
	intDensity[intDensity<0]<-0

min(which(origSequence>=trait),grainScale)


	if (x==bounds[2]){
		X[grainScale] <- 1
	}else{
		nx <- (x-bounds[1])*origIntLength
		ix <- floor(nx)
		ux <- nx-ix
		X[ix+2] <- ux
		X[ix+1] <- 1-ux
		}
	X<-X*origIntLength	
		

		
	


	
	probDivergence <- potentialMatrix %*% X
	# round all up to zero at least
	probDivergence[probDivergence<0] <-	0
	#					
	# sample from this probability distribution
		# to get divergence over a time-step
	newdisplacement <- sample(
		x=origSequence,
		size=1,
		prob= proptemp/sum(proptemp))-states	





	
	
	#subtract current states because we want displacement
    newdisplacement <- rpgm::rpgm.rnorm(n = length(states), mean = (attractor-states)*attraction, sd = sd) 

    return(newdisplacement)
    }	
