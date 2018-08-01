




#' @rdname intrinsicModels
#' @export
landscapeFPK_Intrinsic <- function(params, states, timefrompresent) {
	#a discrete time FPK

		# actualistic bounds on function
	getTraitBounds<-function(x){
		bounds <- c(
			min(x)-((max(x) - min(x))/2),
			min(x)+((max(x) - min(x))/2)
			)
		return(bounds)
		}

	# equation for getting potential	
	potentialFun<-function(x,a,b,c){
		# V(x)=ax4+bx2+cx 
		Vres <- (a*(x^4))+(b*(x^2))+(c*x)
		return(Vres)
		}

	############################################

	# need traits to calculate bounds
	bounds<-getTraitBounds(trait)


	# parameters: a,b,c, rootState, sigma


# lexical scoping example
funA<-function(a,b,c){
	z<-sum(a,b,c)
	newFun<-function(x){
		x*z
		}
	return(newFun)
	}

funB<-funA(a=1,b=2,c=3)
funB(1)
funB(2)
funB<-funA(a=2,b=3,c=4)
funB(1)
funB(2)


	
	# example from vignette for 
		# should make two peak landscape
	params<-c(
		a=10,
		b=0.5,
		c=0,
		sigma=1,
		,bounds)

	a <- params[1]
	b <- params[2]
	c <- params[3]
	sigma <- params[4]
	bounds <- params [5:6]


#	rootState <- 

	# sim controls
	grainscaleFPK = 100
	
	
    sd <- params[1]
    attractor <- params[2]
	#in this model, this should be between zero and one
    attraction <- params[3]    
	
	
	
	
	
	
	
	




##############################################################################

# all of the following only need to be run when the parameters of FPK are changed
	
	# need a function that will return last values returned if all parameters match
	# or will calculate new values if they don't match


# landscape descriptor function
# over the arbitrary interval (-1.5 : 1.5)
seqInterval<-seq(from=-1.5,to=1.5,length.out=grainScale)
# # potentialVector is numeric vector representing the potential
	# length = grainScale
# V(x)=ax4+bx2+cx 
potentialVector<-potentialFun(
	x=seqInterval,
	a=a,b=b,c=c)	
#
# Coefficient of Diffusion of the Model
dCoeff <- log((sigma)^2/2) 
#
# Transition matrix describing prob of evolving between two sites in the trait grid in an infinitesimal time step.	
# Create and diagonalize the transition matrix that has been discretized
# returns: the transition matrix going forward in time, for simulating traits only
#
# make empty matrix
M <- matrix(0,grainScale,grainScale)
#assign values not on outer rows/columns
for (i in 2:(grainScale-1)){
	M[i-1,i] <- exp((potentialVector[i]-potentialVector[i-1])/2)
	M[i+1,i] <- exp((potentialVector[i]-potentialVector[i+1])/2)
	M[i,i] <- -(M[i-1,i]+M[i+1,i])
	}
# assign values to transition matrix not assigned by for loop
M[2,1] <- exp((potentialVector[1]-potentialVector[2])/2)
M[1,1] <- -M[2,1]
M[grainScale-1,grainScale] <- exp((potentialVector[grainScale]-potentialVector[grainScale-1])/2)
M[grainScale,grainScale] <- -M[grainScale-1,grainScale]
# eigen of transition matrix
eig <- eigen(M)
#
passage <- matrix(NA,dim(eig$vectors)[1],dim(eig$vectors)[2])
for (col in 1:dim(eig$vectors)[2]){
	passage[,col] <- Re(eig$vectors[,col])
	}  
# 			
# obtain the matrix diagonal
tPassage <- solve(passage,tol = 1e-30) 
tau <- (((bounds[2]-bounds[1])/(grainScale-1))^2)
vDiag <- diag(Re(eig$values)),
diag_expD <- exp(dCoeff)/tau*diag(vDiag) 


###############################################



###############################################

	x <- states
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
		
	

	# Convolution product over one branch
	expD <- matrix(0,grainScale,grainScale)
	# old time-dep vertsion
	#diag(expD) <- exp(t*diag_expD)
	diag(expD) <- exp(diag_expD)
	a <- passage%*%expD%*%tPassage%*%X
	# prevent rounding errors for small numbers	
	proptemp<-apply(a,1,function(x) max(x,0))
					
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
