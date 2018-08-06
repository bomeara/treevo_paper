

# actualistic bounds on function
getTraitBoundsFPK<-function(x){
	bounds <- c(
		min(x)-((max(x) - min(x))/2),
		max(x)+((max(x) - min(x))/2)
		)
	return(bounds)
	}

getTraitIntervalDensityFPK<-function(trait,origIntLength,
		origSequence,grainScaleFPK){
	# need a vector, length = grainScaleFPK
		# with zeroes in intervals far from current trait value
	# and with density=1 distributed in interval=origIntLength
		# centered around the original trait value
	# for whatever reason, Boucher et al's code
		# chooses the whole interval BEFORE the trait value
		# unclear why you would do that...	
	traitInterval<-c(trait-origIntLength/2, trait+origIntLength/2)
	intDensity<-ifelse(trait>origSequence[-1],
		origSequence[-1]-traitInterval[1],
		traitInterval[2]-origSequence[-grainScaleFPK]
		)
	if(trait==origSequence
	intDensity[intDensity<0]<-0
	if(length(intDensity)!=grainScaleFPK){
		stop("intDensity is not calculated with correct length")
		}
	return(intDensity)
	}




getTraitIntervalDensityFPK<-function(trait,origIntLength,
		origSequence,grainScaleFPK){




	#### example dataset for testing
	# grainScaleFPK<-100
	# origIntLength<-0.23
	# origSequence<-(-20:(-20+grainScaleFPK))*origIntLength
	# trait<-(-1.83)
	# trait<-max(origSequence)
	# trait<-min(origSequence)
	####################################################
	traitRange<-c(trait-origIntLength/2,
		trait+origIntLength/2)
	#
	intDensity<-rep(NA,grainScaleFPK)
	intDensity[2:(grainScaleFPK-1)]<-sapply(2:(grainScaleFPK-1), function(i) max(0,
		min(origSequence[i+1],traitRange[2])-max(origSequence[i],traitRange[1])))
	#
	# special calculations for first and last
	if(traitRange[2]<origSequence[2]){
		intDensity[1]<-origIntLength
	}else{
		intDensity[1]<-0
		}
	#
	if(traitRange[1]>origSequence[grainScaleFPK-1]){
		intDensity[grainScaleFPK]<-origIntLength
	}else{
		intDensity[grainScaleFPK]<-0
		}
	


	#
	# make matrix with interval start and ends
	intBounds<-cbind(
		origSequence[-grainScaleFPK],
		origSequence[-1]
		)
	intDensity<-rep(0,grainScaleFPK)
	for(i in 1:(grainScaleFPK-1)){
		# if the range-start is less than int-end 
			# & range-end is more than int-start
		if(traitRange[1]<intBounds[i,2] 
				& traitRange[2]>intBounds[i,1]){
			# there's overlap, so calculate total overlap
			intDensity[i]<- min(intBounds[i,2],traitRange[2])-max(intBounds[i,1],traitRange[1])
			}
		}


		#
		

		# if the range-start is less than int-end
			# & range-end is greater than int-start - has START IN IT
		if(traitRange[1]>=intBounds[i,1] & traitRange[1]<intBounds[i,2]){
			intDensity[i]<- min(intBounds[i,2],traitRange[2])-traitRange[1]
			}



		}else{
			
			}
		}

	intDensity<-ifelse(trait>origSequence[-1],
		origSequence[-1]-traitRange[1],
		traitRange[2]-origSequence[-grainScaleFPK]
		)
	#if(trait==origSequence
	intDensity[intDensity<0]<-0
	if(length(intDensity)!=grainScaleFPK){
		stop("intDensity is not calculated with correct length")
		}

	#sum(intDensity)==origIntLength

	return(intDensity)
	}




# equation for getting potential under FPK	
potentialFunFPK<-function(x,a,b,c){
	# V(x)=ax4+bx2+cx 
	Vres <- (a*(x^4))+(b*(x^2))+(c*x)
	return(Vres)
	}

#' @rdname 
#' @export
landscapeFPK_Intrinsic <- function(params, states, timefrompresent,	
		grainScaleFPK = 100	# sim controls
		) {
	#
	#a discrete time Fokker–Planck–Kolmogorov model (FPK)
		# V(x)=ax4+bx2+cx 
	# parameters: a,b,c, sigma
		# dependent on former trait value
	#
	# describes a potential surface where height corresponds to
		# tendency to change in that direction
	# allows for multiple optima, different heights to those optima
		#also collapses to BM and OU1
	#
	# From Boucher et al:
	# Finally, note that both BM and the OU model are special cases of the FPK
	# model: BM corresponds to V(x)=0 and OU to
	# V(x)=((alpha/sigma^2)*x^2)-((2*alpha*theta/(sigma^2))*x)
	#
	# following code is loosely based on function Sim_FPK from package BBMV
	#
	# all of the following only need to be run
		# when the parameters of FPK are changed
	# this *could* be pre-calculated for a single run with lexical scoping	
	#	
	# landscape descriptor function
	# over the arbitrary interval (-1.5 : 1.5)
	arbSequence<-seq(from=-1.5,to=1.5,
		length.out=grainScaleFPK)
	#
	# get bounds from params
	bounds <- params[5:6]
	#
	# translate  to original trait scale
	origSequence<-seq(from=bounds[1],to=bounds[2],
		length.out=grainScaleFPK)
	origIntLength<-abs((bounds[2]-bounds[1])/(grainScaleFPK-1))
	# # potentialVector is numeric vector representing the potential
		# length = grainScaleFPK
	# V(x)=ax4+bx2+cx 
	potentialVector<-potentialFunFPK(
		x=arbSequence,
		a=params[1],b=params[2],c=params[3])	
	#
	# Coefficient of Diffusion of the Model
	dCoeff <- log((params[4])^2/2)   # log((sigma^2)/2)
	#
	# Transition matrix describing prob of evolving between two sites in the trait grid in an infinitesimal time step.	
	# Create and diagonalize the transition matrix that has been discretized
	# returns: the transition matrix going forward in time, for simulating traits only
	#
	# make empty matrix
	expD <- tranMatrix <- matrix(0,grainScaleFPK,grainScaleFPK)
	#assign values not on outer rows/columns
	for (i in 1:(grainScaleFPK)){
		if(i>1){
			tranMatrix[i-1,i] <- exp((potentialVector[i]-potentialVector[i-1])/2)
			}
		if(i<grainScaleFPK){
			tranMatrix[i+1,i] <- exp((potentialVector[i]-potentialVector[i+1])/2)
			}
		# rate of staying in place is negative sum of neighboring cells
		neighbors<-c(ifelse(i>1,tranMatrix[i-1,i],0),
				ifelse(i<grainScaleFPK,tranMatrix[i+1,i],0)
				)
		tranMatrix[i,i] <- (-sum(neighbors))
		}
	# eigenvalues and eigenvectors of transition matrix
		# take only real components
	eigTranMatrix <- lapply(eigen(tranMatrix),Re)
	# 			
	# solve the eigenvectors
	solvedEigenvectors <- solve(eigTranMatrix$vectors,tol = 1e-30) 
	#
	# get expected dispersion
	# scale expected dispersion to original trait scale
		# squared distance between points in resolution of trait scale
		# (tau from Boucher et al.'s original code)
	origScaler <- origIntLength^2
	# assign dispersion to diagonal of expD
	diag(expD) <- exp(exp(dCoeff)/origScaler*eigTranMatrix$values)
	# previous time-dep version from Boucher et al's code
		# diag(expD) <- exp(t*diag_expD)
	#
	# take dot product of expD, eigenvectors and solved eigenvectors
	# get matrix of potential for future trait values
		# given potentialVector alone, not yet considering previous values
	potentialMatrix <- eigTranMatrix$vectors%*%expD%*%solvedEigenvectors
	#
	###############################################
	#######################################################
	#
	# need a vector, length = grainScaleFPK
		# with zeroes in intervals far from current trait value
	# and with density=1 distributed in interval=origIntLength
		# centered around the original trait value
	intDensity<-getTraitIntervalDensityFPK(
		trait=states,
		origIntLength=origIntLength,
		origSequence=origSequence,
		grainScaleFPK=grainScaleFPK)
	#
	print(str(intDensity));print(str(potentialMatrix))
	probDivergence <- potentialMatrix %*% intDensity
	# round all up to zero at least
	probDivergence[probDivergence<0] <-	0
	#					
	# sample from this probability distribution
		# to get divergence over a time-step
	newTraitPosition <- sample(
		x=origSequence,
		size=1,
		prob= proptemp/sum(proptemp))
	# subtract the current trait position so to get divergence
	newDisplacement<-newTraitPosition-states
	return(newDisplacement)
	}




traitData<-rnorm(100,0,1)
# need traits to calculate bounds
bounds<-getTraitBoundsFPK(traitData)

# example from vignette for 
	# should make two peak landscape
params<-c(
	a=1,
	b=0.5,
	c=0,
	sigma=1,
	bounds)
trait<-traitData[1]


landscapeFPK_Intrinsic(params=params, states=trait, timefrompresent=NULL)


