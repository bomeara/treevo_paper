
extractIntrinsic_from_prcOut<-function(prcOut){
	res <- list(
		intrinsicFn = prcOut$intrinsicFn, 
		intrinsicValues = prcOut$parMeansList$intrinsic, 
		startingValues = prcOut$parMeansList$starting
		)
	return(res)
	}
	
extractExtrinsic_from_prcOut<-function(prcOut){
	res <- list(
		extrinsicFn = prcOut$extrinsicFn, 
		extrinsicValues = prcOut$parMeansList$extrinsic
		)
	return(res)
	}

runAnalysis <- function(
		runParameters, nSimTrait, ratePriorError,
		#
		anolisTreeList, anolisSize,
		aquilegiaTreeList, aquilegiaSpurLength,
		idealTrees,
		#
		indepAnalyses_intrinsicOut, 
		indepAnalyses_extrinsicOut,
		#
		# presets
		generation.time,
		multicore,
		coreLimit,				
		nRuns,
		nStepsPRC,
		numParticles, 
		nInitialSimsPerParam,
		StartSims	
		){
	################
	################################################
	# define MCMC / ABC control parameter list
	controlsList <- list(
		# standard controls, don't need to be changed
		standardDevFactor=0.2,				
		epsilonProportion=0.7,
		epsilonMultiplier=0.7,
		stopRule = FALSE, 
		plot = FALSE,
		verboseParticles = FALSE,	
		#
		# controls that may need to be changed
		generation.time = generation.time,
		multicore = multicore,
		coreLimit = coreLimit,				
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
	if(runParameters$prior != "rexp_with_mean_NOT_at_true_sigmasq"){
		# then do *NOT* apply the error to the sigmasq prior
		# Reset ratePriorError to 1
		ratePriorError <- 1
		}
	##########################################
	#
	#
	# doRun.Intrinsic
	# 
	if(runParameters$doRun.Intrinsic == "Pure_BM"){
		intrinsicFunctionToFit <- brownianIntrinsic
		#
		intrinsicArgList <- list(
			intrinsicPriorsFns = c("exponential"), 
			intrinsicPriorsValues = list(10 * ratePriorError)
			)
		}
	#
	if(runParameters$doRun.Intrinsic == "BM_LowerBound"){
		intrinsicFunctionToFit <- boundaryMinIntrinsic
		#
		intrinsicArgList <- list(
			intrinsicPriorsFns=c("exponential", "normal"),
			intrinsicPriorsValues=list(10, c(-10, 1))
			)
		}
	#
	if(runParameters$doRun.Intrinsic == "3Opt2Bound"){
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
				"normal", "normal", "normal"
				),
			intrinsicPriorsValues=list(
				10, 10, 1, 
				c(15,20),c(20,25),
				c(10,1), c(20,1), c(30,1)
				)
			)
		}
	#
	if(runParameters$doRun.Intrinsic == "Time_AutoRegressive_Model"){
		intrinsicFunctionToFit <- autoregressiveIntrinsic
		#
		intrinsicArgList <- list(
			intrinsicPriorsFns=c("exponential", "normal"),
			intrinsicPriorsValues=list(10, c(-10, 1))
			)
		}
	#
	# doRun.Extrinsic
	#
	if(runParameters$doRun.Extrinsic =="Null"){
		extrinsicFunctionToFit <- nullExtrinsic
		#
		extrinsicArgList <- list(
			extrinsicPriorsFns = c("fixed"), 
			extrinsicPriorsValues = list(0)
			)
		}
	#
	if(runParameters$doRun.Extrinsic == "Displacement"){
		extrinsicFunctionToFit  <- ExponentiallyDecayingPushExtrinsic
		#
		extrinsicArgList <- list(
			extrinsicPriorsFns = c("exponential", "normal", "exponential"), 
			# \code{ExponentiallyDecayingPushExtrinsic}
			# with parameters \code{params = sd, maximumForce, halfDistance}
			extrinsicPriorsValues = list(10, c(1,1), 10)
			)
		}
	#########################################
	# nTraitSetsPerSimTree is 1 unless empiricalTraitData is
		# "SIMULATED" in which case it is nSimTrait
	nTraitSetsPerSimTree <- 1
	if(runParameters$empiricalTraitData == "SIMULATED"){
		nTraitSetsPerSimTree <- nSimTrait
		}
	#
	#########################################
	#
	#treeSet
	#
	# if the treeSet is "Ideal-Simulated"
	# then the number of simulated tree types and 
		# number of tip-totals per simulated tree type is 3, other 1
	nSimTreeTypes <- nTipNumbersPerSimTreeType <- 1
	#
	if(runParameters$treeSet == "empirical_anolis_tree"){
		treeList <- anolisTreeList
		}
	#
	if(runParameters$treeSet == "empirical_Aquilegia_tree"){
		treeList <- aquilegiaTreeList
		}
	#
	if(runParameters$treeSet == "Ideal_Simulated"){
		treeList <- idealTrees
		nSimTreeTypes <- nTipNumbersPerSimTreeType <- 3
		}
	#
	################################################	
	# need to make trait data for every tree in treeList	
	#
	# traitDataList will be a list with each element corresponding to a tree
	# and sub list corresponding to trait data to be analyzed on that tree
	#
	traitDataList <- list()
	#
	for (tree_i in 1:length(treeList)){
		#
		# empiricalTraitData
		#
		if(runParameters$empiricalTraitData == "Anolis_Size_Data"){
			# need a list of trait sets (of length 1)
			traitDataList[[tree_i]]  <- anolisSize
			}
		#
		if(runParameters$empiricalTraitData == "Aquilegia_Nectar_Spur_Data"){
			# need a list of trait sets (of length 1)
			traitDataList[[tree_i]]  <- aquilegiaSpurLength
			}
		#
		if(runParameters$empiricalTraitData == "SIMULATED"){
			#
			simTrait.Intrinsic <- runParameters$simTrait.Intrinsic
			simTrait.Extrinsic <- runParameters$simTrait.Extrinsic
			#
			# simTrait.Intrinsic
			# ALSO need estimates of parameters from previous analyses needed for later simulations
			# 	 need to make part of output from doRun if not already
			#
			if(is.na(simTrait.Intrinsic)){
				stop("The intrinsic model for a simulated trait dataset is given as NA")
			}else{
				# call respective analysis, take parameters from it
				simTraitIntrinsicArgs <- list(
					intfn = indepAnalyses_intrinsicOut[[simTrait.Intrinsic]]$intrinsicFn,
					intPar = indepAnalyses_intrinsicOut[[simTrait.Intrinsic]]$intrinsicValues,
					startPar = indepAnalyses_intrinsicOut[[simTrait.Intrinsic]]$startingValues
					)			
				}
			#
			# simTrait.Extrinsic
			#
			if(is.na(runParameters$simTrait.Extrinsic)){
				stop("The extrinsic model for a simulated trait dataset is given as NA")
			}else{
				if(simTrait.Extrinsic == "Null"){
					simTraitExtrinsicArgs <- list(
						extfn = nullExtrinsic,
						extPar = c(0), 
						)
				}else{
					simTraitExtrinsicArgs <- list(
						extfn = indepAnalyses_extrinsicOut[[simTrait.Extrinsic]]$extrinsicFn,
						extPar = indepAnalyses_extrinsicOut[[simTrait.Extrinsic]]$extrinsicValues
						)	
					}		
				}	
			#####################
			# now have to simulate traits
				# save to the list of trait sets
			traitDataList[[tree_i]]  <- doSimulation(
				phy = treeList[[tree_i]], 
				intrinsicFn = simTraitIntrinsicArgs$intFn, 
				extrinsicFn = simTraitExtrinsicArgs$exFn, 
				startingValues = simTraitIntrinsicArgs$startPar,
				intrinsicValues = simTraitIntrinsicArgs$intPar, 
				extrinsicValues = simTraitExtrinsicArgs$exPar, 
				generation.time = generation.time
				)	
			}
		}
	#################################################
	#
	# nDoRun
	#
	# calculate the number of doRun statements for this analysis-run
	# product of nSimTreeTypes and nTipNumbersPerSimTreeType and nSimTrait
	nDoRun <- nSimTreeTypes * nTipNumbersPerSimTreeType * nTraitSetsPerSimTree
	# should be one 1, 10 or 90... probably
	#	
	##########################################################
	# now run doRun across trees, trait datasets
	#
	for (trait_j in length(traitDataList)){
		# define job name
		jobNameRun <- paste0(
			runParameters$runLabel,
			"_", format(Sys.time(), "%m-%d-%y"),
			"_", trait_j
			)
		#
		traitDataToUseForThisRun <- traitDataList[[trait_j]]
		#if(is.list(traitDataToUseForThisRun)){
		#	traitDataToUseForThisRun <- unlist(traitDataToUseForThisRun)
		#	}
		#print(traitDataToUseForThisRun)
		#x<-getBM(treeList[[trait_j]], 
		#	trait = traitDataToUseForThisRun)
		#
		argListForDoRun <- list(
			##############
			phy = treeList[[trait_j]],
			traits = traitDataToUseForThisRun,
			#
			intrinsicFn = intrinsicFunctionToFit,
			extrinsicFn = extrinsicFunctionToFit,
			#
			startingPriorsFns = "normal", 
			startingPriorsValues = 
				list(c(
					mean(traitDataToUseForThisRun),
					sd(traitDataToUseForThisRun)
					)), 
			#########
			#
			intrinsicPriorsFns =
				intrinsicArgList$intrinsicPriorsFns, 
			intrinsicPriorsValues =
				intrinsicArgList$intrinsicPriorsValues, 
			#########
			#
			extrinsicPriorsFns = 
				extrinsicArgList$extrinsicPriorsFns, 
			extrinsicPriorsValues = 
				extrinsicArgList$extrinsicPriorsValues, 
			#########
			#
			jobName = jobNameRun
			)
		#
		# add control arguments
		argListForDoRun <- c(argListForDoRun, controlsList)
		#	
		doRun_out <- do.call(what = doRun_prc,
			# arguments
			args = argListForDoRun
			)
		}
	#
	return(doRun_out)
	}
	