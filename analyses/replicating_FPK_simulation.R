

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
	a=10,
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
seqInterval<-seq(from=-1.5,to=1.5,length.out=grainScale)
# # potentialVector is numeric vector representing the potential
	# length = grainScale
# V(x)=ax4+bx2+cx 
potentialVector<-potentialFun(
	x=seqInterval,
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


# Convolution product over one branch
# what is this tau? depends on bounds, grainScale
tau <- (((bounds[2]-bounds[1])/(grainScale-1))^2)
# old time-dep version from Boucher et al's code
#diag(expD) <- exp(t*diag_expD)
diag(expD) <- exp(exp(dCoeff)/tau*eigM$values)


###############################################



###############################################

	x <- trait
	#t <- tree$edge.length[i]

	# write to which point of the grid a given position belongs to, 'continuous' version
	X <- rep(0,grainScale)  
	
	


	if (x==bounds[2]){
		X[grainScale] <- 1
	}else{
		nx <- (grainScale-1)*(x-bounds[1])/(bounds[2]-bounds[1])
		ix <- floor(nx)
		ux <- nx-ix
		X[ix+2] <- ux
		X[ix+1] <- 1-ux
		}

	X<-X*(grainScale-1)/(bounds[2]-bounds[1])	
		
	# get the probability density / likelihood of each future trait value
		# given potentialVector and current trait value
	# propagate the trait forward in time
		
	


	proptemp2 <- passage%*%expD%*%solvedEigenvectors
	proptemp2%*%X
	proptemp <- passage%*%expD%*%solvedEigenvectors%*%X
	proptemp

	# prevent rounding errors for small numbers	
	proptemp<-apply(proptemp,1,function(x) max(x,0))
					
	# sample from this probability distribution	to get a descendent node value
	newdisplacement <- sample(
		x=seq(
			from=bounds[1],
			to=bounds[2],
			length.out=grainScale
			)-,
		size=1,
		prob= proptemp/sum(proptemp))-states	





	
	
	#subtract current states because we want displacement
    newdisplacement <- rpgm::rpgm.rnorm(n = length(states), mean = (attractor-states)*attraction, sd = sd) 

    return(newdisplacement)
    }	
