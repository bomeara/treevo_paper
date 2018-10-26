###################################################
# Individual Empirical Analyses and Simulations
##################################################

# Control Box



# number of simulated trait datasets to do for simulated-trait runs
nSimTrait <- 10

# error for run with mis-specified prior on sigmasq in a pure-BM model
	# the mean of the normal prior is multiplied by this value
	# 100 = mean of rate prior is off by two orders of magnitude!
ratePriorError <- 100

# control parameters for multicore and simulation resolution
multicore <- TRUE 
coreLimit <- 6
generation.time <- 10000 

# control parameters for MCMC / ABC
nRuns <- 2 
nStepsPRC <- 3 
numParticles <- 20 
nInitialSimsPerParam <- 10 
StartSims <- 10

##################################################

setwd("d://dave//workspace//treevo_paper//")

library(ape)
library(TreEvo)

######################################
# get empirical data
	
# 1) Anolis
	# repulsion - adaptive landscape dynamics - multi optima	

# obtain anolis tree - from Poe et al. 2017 (SystBiol)
	# their time-tree
anolisTree<-read.tree(file="datasets//anolis_PoeEtAl2018_datedMCCw0.5burnin.tre")
# make into a multiPhylo list
anolisTreeList <- list(anolisTree = anolisTree)
class(anolisTreeList) <- "multiPhylo"

# obtain anolis trait data - 
	# Snout-Vent body-size data from Poe et al. 2018 (AmNat)
anolisTrait<-read.table("datasets//anolis_lntraits_matched_tabdelim_07-24-18.txt",
	header=TRUE,row.names=1)
anolisSize<-anolisTrait[,1]

# 2) Aquilegia
	# whittall et al. model of nectar spur increase in size

# obtain aquilegia tree (from Whittall and Hodges 2007?)
aquilegiaTree<-read.tree("datasets//aquilegia_Whttall&Hodges2007_figuredMCC.tre")
# make into a multiPhylo list
aquilegiaTreeList <- list(aquilegiaTree = aquilegiaTree)
class(aquilegiaTreeList) <- "multiPhylo"

# obtain aquilegia trait data (from Whittall and Hodges 2007?)
	# need both nectur spur lengths and regime data
	#
aquilegiaTrait<-read.table("aquilegia_traitData.txt", header=FALSE, row.names=1)

# get just nectur spur length
aquilegiaSpurLength<-aquilegiaTrait[,2]
# and take the natural log
	# (note that the third column of the table was already the natural log)
	# previous code from Brian had 'log(data[,3])' - log of a log
aquilegiaSpurLength<-log(aquilegiaSpurLength)

# legacy aquilegia code from Brian O'Meara:
# 
# assume generation time of 10 years (its a perennial plant), 
	# following Cooper et al. Plos ONe 2010 
	# Genetic Variation at Nuclear loci fails to distinguish group is about 3 MY,
		# phy height is 3. So each unit = 1,000,000 years or thus 100,000 generations
# TreeYears=100000
# timeStep<-1/TreeYears
# totalTreeLength=TreeYears*sum(phy$edge.length) #how many generations are represented
# number of expected polinator shifts based on parsimony is 7:
# parsimonyShifts=7
# pollinatorShiftRate=parsimonyShifts/totalTreeLength

# aquilegia regimes - pollinator syndromes
aquilegiaPollinators<-aquilegiaTrait[,14]
# regimes coded 0, 1, 2
	# 0 is bumble-bee, 1 is humming-bird, 2 is hawkmoth
# this probably won't be used directly?
# could use for post-analysis comparisons? Hmm
	

###############################################################################
# generate sets of ideal trees for doing simulations on
   # idealTreeSets = c("Ideal-Balanced", "Ideal-Pectinate", "Ideal-Star") 
   # nTipSets = c(8, 16, 64)   
   
idealTrees<-list(
    #
    balanced_n8 = stree(n=8, type = "balanced", tip.label = NULL),
    balanced_n16 = stree(n=16, type = "balanced", tip.label = NULL),
    balanced_n64 = stree(n=64, type = "balanced", tip.label = NULL),
    #
    pectinate_n8 = stree(n=8, type = "left", tip.label = NULL),
    pectinate_n16 = stree(n=16, type = "left", tip.label = NULL),
    pectinate_n64 = stree(n=64, type = "left", tip.label = NULL), 
    #
    star_n8 = stree(n=8, type = "star", tip.label = NULL),
    star_n16 = stree(n=16, type = "star", tip.label = NULL),
    star_n64 = stree(n=64, type = "star", tip.label = NULL)
    )
# make multiPhylo	
class(ideaTrees)<-"multiPhylo"
# compress tip labels? No, I don't think that works for trees of different sizes
	#	trees<-.compressTipLabel(trees)
######################################################################################
# time to get table, process the inputs listed
# 		
# get simulation run table
simRunTable<-read.csv(header=TRUE,
	stringsAsFactors=FALSE,
	file="simulation_sets_parameters_table.csv")	
#
# number of analyses
nAnalyses<-nrow(simRunTable)	
#
# which analyses are independent or dependent
whichIndependentPrevRun<-which(!as.logical(simRunTable$dependentPrevRun))	
whichDependentPrevRun<-which(as.logical(simRunTable$dependentPrevRun))
#
# create list for saving analysis output 
analysisOutput <- list()
names(analysisOutput)<-1:nAnalyses
#
# Let's runs the analyses! 
#
# run all independent analyses
for (i in whichIndependentPrevRun){
	analysisOutput[[i]]<-
	}
#






# run all dependent analyses
for (i in whichDependentPrevRun){
	analysisOutput[[i]]<-
	}	
		
		
		
		
# inputs needed from script	above
nSimTrait
ratePriorError

anolisTreeList
aquilegiaTreeList

anolisSize
aquilegiaSpurLength

idealTrees


	generation.time=generation.time,
	multicore=multicore,
	coreLimit=coreLimit,				
	nRuns = nRuns, 
	nStepsPRC = nStepsPRC, 
	numParticles = numParticles, 
	nInitialSimsPerParam = nInitialSimsPerParam, 
	StartSims = StartSims	


		
	generation.time
	multicore
	coreLimit			
	nRuns  
	nStepsPRC 
	numParticles 
	nInitialSimsPerParam 
	StartSims 	





		
		

################################################
# define MCMC / ABC control parameter list
controlsList <- list(
	# standard controls, don't need to be change
	standardDevFactor=0.2,				
	epsilonProportion=0.7,
	epsilonMultiplier=0.7,
	stopRule = FALSE, 
	plot=FALSE,
	verboseParticles=FALSE,	
	#
	# controls that may need to be changed
	generation.time=generation.time,
	multicore=multicore,
	coreLimit=coreLimit,				
	nRuns = nRuns, 
	nStepsPRC = nStepsPRC, 
	numParticles = numParticles, 
	nInitialSimsPerParam = nInitialSimsPerParam, 
	StartSims = StartSims	
	)
#
##################################################	
# rate prior error
# sigmasq state prior has an error if "rexp_with_mean_NOT_at_true_sigmasq"
# 
if (prior != "rexp_with_mean_NOT_at_true_sigmasq"){
	# then do *NOT* apply the error to the sigmasq prior
	ratePriorError <- 1
	}
##########################################
#
#
# doRun.Intrinsic
# 
if (doRun.Intrinsic == "Pure_BM"){
	intrinsicFunctionToFit <- brownianIntrinsic
	#
	intrinsicArgList <- list(
		intrinsicPriorsFns = c("exponential"), 
		intrinsicPriorsValues = list(10 * ratePriorError)
		)
	}
#
if (doRun.Intrinsic == "BM_LowerBound"){
	intrinsicFunctionToFit <- boundaryMinIntrinsic
	#
	intrinsicArgList <- list(
		intrinsicPriorsFns=c("exponential","normal"),
		intrinsicPriorsValues=list(10, c(-10, 1))
		)
	}
#
if (doRun.Intrinsic == "3Opt2Bound"){
	intrinsicFunctionToFit <- multiOptima3IntrinsicMaxBoundary2
	#
	intrinsicArgList <- list(
		# breakdown of params:
			# params[1] is dispersion (sigma)
			# params[2] is alpha (strength of attraction to an optima)
			# params[3] is rho, an exponent scaling the weighting of distance to optima
				# this parameter will control switching optima
			# params[4:5] is the max boundary, for the two lower regimes regimes
			# params[6:8] describes theta (optima) values for each of the three regimes
		intrinsicPriorsFns=c(
			# we'll make rate an exponential prior, rate 10
			"exponential",
			# we'll make alpha an exponential prior, rate 10
			"exponential",
			# well make rho an exponential, rate 1
			"exponential",
			# let's place bounds based on whitall and hodges:
			"uniform","uniform",
			# normal priors for optima
			"normal","normal","normal"
			),
		intrinsicPriorsValues=list(
			10, 10, 1, 
			c(15,20),c(20,25),
			c(10,1), c(20,1), c(30,1)
			)
		)
	}
#
if (doRun.Intrinsic == "Time_AutoRegressive_Model"){
	intrinsicFunctionToFit <- autoregressiveIntrinsic
	#
	intrinsicArgList <- list(
		
		intrinsicPriorsFns=c("exponential","normal"),
		intrinsicPriorsValues=list(10, c(-10, 1))
		)
	}
#
# doRun.Extrinsic
#
if (doRun.Extrinsic =="Null"){
	extrinsicFunctionToFit <- nullExtrinsic
	#
	extrinsicArgList <- list(
		extrinsicPriorsFns = c("fixed"), 
		extrinsicPriorsValues = list(0)
		)
	}
#
if (doRun.Extrinsic =="Displacement"){
	extrinsicFunctionToFit <-ExponentiallyDecayingPushExtrinsic
	#
	extrinsicArgList <- list(
		extrinsicPriorsFns = c("exponential","normal","exponential"), 
		# \code{ExponentiallyDecayingPushExtrinsic} with parameters \code{params = sd, maximumForce, halfDistance}
		extrinsicPriorsValues = list(10,c(1,1),10)
		)
	}
#########################################
# nTraitSetsPerSimTree is 1 unless empiricalTraitData is "SIMULATED" in which case it is nSimTrait
nTraitSetsPerSimTree<-1
if(empiricalTraitData == "SIMULATED"){
	nTraitSetsPerSimTree<-nSimTrait
	}
#
#########################################
#
#treeSet
#
# if the treeSet is "Ideal-Simulated"
# then the number of simulated tree types and 
	# number of tip-totals per simulated tree type is 3, other 1
nSimTreeTypes<-nTipNumbersPerSimTreeType<-1
#
if(treeSet == "empirical_anolis_tree"){
	treeList <- anolisTreeList
	}
#
if(treeSet == "empirical_Aquilegia_tree"){
	treeList <- aquilegiaTreeList
	}
#
if(treeSet=="Ideal_Simulated"){
	treeList<-idealTrees
	nSimTreeTypes<-nTipNumbersPerSimTreeType<-3
	}
####################################
#
# nDoRun
#
# calculate the number of doRun statements for this analysis-run
# product of treeTypes and nTipNumbersPerSimTreeType and nSimTrait
nDoRun <- nSimTreeTypes * nTipNumbersPerSimTreeType * nTraitSetsPerSimTree
# should be one 1, 10 or 90... probably
#	
################################################	
# need to make trait data for every tree in treeList	
#




for (tree_i in 1:length(treeList)){
#
# traitDataList will be a list with each element corresponding to a tree
# and sub list corresponding to trait data to be analyzed on that tree
#
traitDataList<-list()
# empiricalTraitData
#
if(empiricalTraitData == "Anolis_Size_Data"){
	# need a list of trait sets (of length 1)
	traitDataList[[tree_i]] <-list(anolisSize = anolisSize)
	}
#
if(empiricalTraitData == "Aquilegia_Nectar_Spur_Data"){
	# need a list of trait sets (of length 1)
	traitDataList[[tree_i]] <-list(aquilegiaSpurLength = aquilegiaSpurLength)
	}
#
if(empiricalTraitData == "SIMULATED"){
	#
	# simTrait.Intrinsic
	# ALSO need estimates of parameters from previous analyses needed for later simulations
	# 	 need to make part of output from doRun if not already
	#
	if(is.na(simTrait.Intrinsic)){
		stop("The intrinsic model for a simulated trait dataset is given as NA")
	}else{
		# ANOLIS BASED MODELS
		if(simTrait.Intrinsic == "An_Emp_BrownMotion"){
			simTraitIntrinsicArgs <- list(
				intfn = brownianIntrinsic,
				intPar = An_Emp_BrownMotion$parMeansList$intrinsic, #whatever run is An_Emp_BrownMotion
				startPar = An_Emp_BrownMotion$parMeansList$starting
				)
			
				

	 
	 

# need function that simply outputs a list with those parameters
	# (median?) expectations from the last MCMC generation
# this function would be run with doRun such that doRun would include with output
# these parameter estimates would be given as a list
	# split into 3 vectors: starting/intrinsic/extrinsic parameters
	# formatted for immediate use as parameter estimates for doSimulation
	# with matching intrinsic/extrinsic functions
	

	 
				
				)
			}	
		#
		if(simTrait.Intrinsic == "An_Emp_Disp"){
			
			}	
		#
		if(simTrait.Intrinsic == "An_Emp_DispBound"){
			
			}	
		#
		if(simTrait.Intrinsic == "An_Emp_Bound"){
			
			}	
		#
		if(simTrait.Intrinsic == "An_Emp_Bound_BoundByStartingState"){
			
			}	
		#
		if(simTrait.Intrinsic == "An_Emp_Bound_BoundByMinValue"){
			
			}	
		#
		if(simTrait.Intrinsic == "An_Emp_Bound_BoundOneRangeAway"){
			
			}	
		#
		if(simTrait.Intrinsic == "An_Emp_TimeReg"){
			
			}	
		#
		if(simTrait.Intrinsic == "Aq_Emp_3Opt2Bound"){
			
			}	
		#
		if(simTrait.Intrinsic == "Aq_Emp_BrownMotion"){
			
			}	
		#
				
			
		}
	#
	# simTrait.Extrinsic
	#
	if(is.na(simTrait.Extrinsic)){
		stop("The extrinsic model for a simulated trait dataset is given as NA")
	}else{
		
		if(simTrait.Extrinsic == "Null"){
			simTraitExtrinsicArgs <- list(
				extfn = nullExtrinsic,
				extPar = c(0), 
				)
			}
		#
		if(simTrait.Extrinsic == "An_Emp_Disp"){
			simTraitExtrinsicArgs <- list(
				extfn = ExponentiallyDecayingPushExtrinsic,
				extPar = anolisBMrun$parMeansList$extrinsic, #whatever run is An_Emp_BrownMotion 
				
				
				
				)			
			}
		#
		if(simTrait.Extrinsic == "An_Emp_Disp"){
			
			}
		#
		if(simTrait.Extrinsic == "An_Emp_DispBound"){
			
			}
		#
	
		
		}
	#####################
	# now have to simulate traits
	simChar <- doSimulation(
		phy = treeList[[tree_i]], 
		
		intrinsicFn = simTraitIntrinsicArgs$intFn, 
		extrinsicFn = exFn, 
		startingValues = simTraitIntrinsicArgs$startPar, #root state
		intrinsicValues = simTraitIntrinsicArgs$intPar, 
		extrinsicValues = c(0), 
		generation.time = generation.time
		)	
	
	# save as a list of trait sets
	traitDataList <- 
	}

# now run doRun across trees, trait datasets
#

	
	

	for (trait_j in length(traitDataList)){
	
		# define job name
		jobNameRun <- 
		#
		traitDataToUse <- traitDataList [[trait_j]]
		#	
		doRun_out <- do.call(what = doRun_prc,
			# arguments
			args = c(
				phy = ,
				traits = traitDataToUse,
				intrinsicFn =
				extrinsicFn = 
				
				intrinsicArgList,
				extrinsicArgList,
				
				# starting state prior
				startingPriorsFns = "normal",
				startingPriorsValues = matrix(
					c(mean(traitDataToUse[, 1]) , sd(traitDataToUse[, 1]))),
				#
				#
				jobName = paste(),
				
				# give control arguments
				controlsList
				)
			)
		}	
	}	

		


  
  

  
  


  )	
		
		
		
		
	
	
###################################################
###### TWO tests of treevo
#################################################

# Unresolved question: Number of particles? Number of generations?

#########################
### (1) test basic BM
	# fixed sigma square (rate)

	#Two different BM priors: unif with truth at say 25th percentile, rexp with mean not at true value (just to be realistic). Similar for root state.

	
#############################	
### (2) interaction model - repulsion between species with max bound 
	# assume log transformed traits, so no min bound
	# use realistic tree set

# Some models with 
	# no actual repulsion
	# moderate (species cross in trait space) 
	# high (no touching happens). 
# Models vary distance of max to the trait data
	# max is very close (most species bounce off it based on starting values), 
	# max is moderately far (start hitting in last 25% of sim), 
	# very far (never hit)
# That’s nine different parameter values, 
	# all of which I guess are scaled effective to the BM rate. 
	# Could do it as six: 
		#1) moderate repulsion, max bound is very close
		#2) moderate repulsion, max bound is moderately close
		#3) moderate repulsion, max bound is very far
		#4) no repulsion, max bound is moderately close
		#5) high repulsion, max bound is moderately close 		
		
#########################
# 3) Time - autoregressive model with optimum based on a factor that changes through time (like O2 concentration)
 # empirical tree, like in simulations

 # Parameters of the regression: 
	# a single variable function to convert O2 to optimal gill size or whatever
	# strength of pull
	# BM wiggle
	
# abd presumably the tracked env factor will be analyzed many times over





	
	

	
	
