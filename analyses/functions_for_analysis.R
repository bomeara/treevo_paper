
extractIntrinsic_from_prcOut<-function(prcOut){
	res <- list(
		intrinsicFn = prcOut[[1]]$intrinsicFn, 
		intrinsicValues = prcOut[[1]]$parMeansList$intrinsic, 
		startingValues = prcOut[[1]]$parMeansList$starting
		)
	return(res)
	}
	
extractExtrinsic_from_prcOut<-function(prcOut){
	res <- list(
		extrinsicFn = prcOut[[1]]$extrinsicFn, 
		extrinsicValues = prcOut[[1]]$parMeansList$extrinsic
		)
	return(res)
	}
	
cleanSimTraitData <- function(simulatedTraitData){
	res <- simulatedTraitData$states
	names(res) <- rownames(simulatedTraitData)
	return(res)
	}				
					

runAnalysis <- function(
		runParameters, 
		nSimTrait, 
		ratePriorError,
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
		nInitialSims,	
		saveData,
		verboseParticles
		){
	##################################################	
	# rate prior error
	# sigmasq state prior has an error
		# if "rexp_with_mean_NOT_at_true_sigmasq"
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
	#
	# treeSet and nDoRun
	#
	# calculate the number of doRun statements for each analysis-run
		# should be one 1, 10 or 90... probably
		# 1 if empirical tree, empirical trait data
		# 10 if empirical tree and 10 simulated trait datasets
		# 90 if 9 simulated trees and 10 simulated trait datasets for *each* tree
	#
	#
	# if the treeSet is "Ideal-Simulated"
		# then the number of simulated tree types and 
		# number of tip-totals per simulated tree type is 3, other 1
	# number of analyses also depends on if empiricalTraitData is "SIMULATED"
	
	#
	if(runParameters$treeSet == "empirical_anolis_tree"){
		treeList <- anolisTreeList
		#		
		if(runParameters$empiricalTraitData == "Anolis_Size_Data"){
			# need a two-level list of trait sets (one tree, one trait dataset)
			traitDataList  <- list(list(anolisSize = anolisSize))
			#
			message(paste0(
				"Performing  an single analysis with the empirical Anolis phylogeny,\n",
				"   and empirical Anolis size trait data."
				))		
			}
		if(runParameters$empiricalTraitData == "SIMULATED"){
			# nDoRun equal to nSimTrait
			nDoRun <- nSimTrait
			#
			message(paste0(
				"Performing ", nDoRun, 
				" analyses on the empirical Anolis phylogeny, with ", 
				nDoRun, 
				" seperately-simulated trait datasets."
				))				
			}			
		}
	#
	if(runParameters$treeSet == "empirical_Aquilegia_tree"){
		treeList <- aquilegiaTreeList
		#		
		if(runParameters$empiricalTraitData == "Aquilegia_Nectar_Spur_Data"){
			# need a two-level list of trait sets (one tree, one trait dataset)
			traitDataList  <- list(list(aquilegiaSpurLength = aquilegiaSpurLength))
			#
			message(paste0(
				"Performing  an single analysis with the empirical Aquilegia phylogeny,\n",
				"   and empirical Aquilegia Nectar Spur trait data."
				))		
			}
		if(runParameters$empiricalTraitData == "SIMULATED"){
			# nDoRun equal to nSimTrait
			nDoRun <- nSimTrait
			#
			message(paste0(
				"Performing ", nDoRun, 
				" analyses on the empirical Aquilegia phylogeny, with ", 
				nDoRun, 
				" seperately-simulated trait datasets."
				))				
			}	
		}
	#
	if(runParameters$treeSet == "Ideal_Simulated"){
		treeList <- idealTrees
		#
		# calculate the number of doRun statements for this analysis-run
		# product of nSimTreeTypes (3) and nTipNumbersPerSimTreeType (3) and nSimTrait
		nDoRun <- nSimTrait * 3 * 3		
		message(paste0(
			"Performing ", nDoRun, 
			" analyses for ", 3,
			" simulated 'idealized' phylogeny classes, with ", 3, 
			" sets of tip values each,\n",
			"    and ", nSimTrait, 
			" simulate trait datasets for each tree."
			))		
		#
		}
	#
	#####################################################################
	# Get Simulated Trait Data
	#
	# traitDataList will be a list with each element corresponding to a tree
		# and sub list corresponding to trait data to be analyzed on that tree
	#
	if(runParameters$empiricalTraitData == "SIMULATED"){
		# go through each tree from tree list
		# iterate, generate nSimTrait simulated trait datasets
		# produce a list where each item represents a tree
			# and each subitem of each list item is a trait dataset
		#
		# first make the empty list
		#
		traitDataList <- list()		
		#
		simTrait.Intrinsic <- runParameters$simTrait.Intrinsic
		simTrait.Extrinsic <- runParameters$simTrait.Extrinsic
		#
		# simTrait.Intrinsic
		# ALSO need estimates of parameters from
		#    previous analyses needed for later simulations
		# 	 need to make part of output from doRun if not already
		#
		if(is.na(simTrait.Intrinsic)){
			stop("The intrinsic model for a simulated trait dataset is given as NA")
		}else{
			# call respective analysis, take parameters from it
			simTraitIntrinsicArgs <- list(
				intrinsicFn = indepAnalyses_intrinsicOut
					[[simTrait.Intrinsic]]$intrinsicFn,
				intrinsicValues = indepAnalyses_intrinsicOut
					[[simTrait.Intrinsic]]$intrinsicValues,
				startingValues = indepAnalyses_intrinsicOut
					[[simTrait.Intrinsic]]$startingValues
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
					extrinsicFn = nullExtrinsic,
					extrinsicValues = c(0)
					)
			}else{
				simTraitExtrinsicArgs <- list(
					extrinsicFn = indepAnalyses_extrinsicOut
						[[simTrait.Extrinsic]]$extrinsicFn,
					extrinsicValues = indepAnalyses_extrinsicOut
						[[simTrait.Extrinsic]]$extrinsicValues
					)	
				}		
			}	
		#####################
		# now have to simulate traits
			# save to the list of trait sets
		#
		for(tree_i in 1:length(treeList)){
			#
			traitDataThisTree <- list()
			#
			for(trait_i in 1:nSimTrait){
				#
				simulatedTraitData  <- doSimulation(
					phy = treeList[[tree_i]], 
					intrinsicFn = simTraitIntrinsicArgs$intrinsicFn, 
					extrinsicFn = simTraitExtrinsicArgs$extrinsicFn, 
					startingValues = simTraitIntrinsicArgs$startingValues,
					intrinsicValues = simTraitIntrinsicArgs$intrinsicValues, 
					extrinsicValues = simTraitExtrinsicArgs$extrinsicValues, 
					generation.time = generation.time
					)	
				#
				traitDataThisTree[[trait_i]] <- cleanSimTraitData(simulatedTraitData)
				}
			#
			traitDataList[[tree_i]]  <- traitDataThisTree
			}
		}
	#############################################################
	##########################################################
	# now run doRun across each trees and its trait datasets
	#
	# first make the empty list for output - two levels!	
	doRun_out <- list()
	#
	for (i in 1:length(treeList)){
		# first iterate over trees
		#
		treeToUse <- treeList[[i]]
		#	
		# empty list
		doRun_out_ThisTree <- list()
		#
		for (j in 1:length(traitDataList[[i]])){
			#
			# define job name
			jobNameRun <- paste0(
				runParameters$runLabel,
				"_tree_",i,
				"_trait_",j,
				"_", format(Sys.time(), "%m-%d-%y")
				)
			#
			traitDataToUseForThisRun <- traitDataList[[i]][[j]]
			#
			#print(traitDataToUseForThisRun)
			#
			doRun_out_ThisTree[[j]] <- doRun_prc(
				##############
				phy = treeToUse,
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
				###########################
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
				jobName = jobNameRun,
				#
				################################################
				# define MCMC / ABC control parameter list
				#
				# controls that may need to be changed
				generation.time = generation.time,
				multicore = multicore,
				coreLimit = coreLimit,				
				#
				numParticles = numParticles, 
				nStepsPRC = nStepsPRC, 
				nRuns = nRuns, 
				nInitialSims = nInitialSims,
				nInitialSimsPerParam = nInitialSimsPerParam, 
				#
				saveData = saveData, 
				verboseParticles = verboseParticles,
				#
				#
				# standard controls, don't need to be changed
				standardDevFactor = 0.20, 
				epsilonProportion = 0.7, 
				epsilonMultiplier = 0.7, 
				#
				validation = "CV", 
				scale = TRUE, 
				variance.cutoff = 95, 
				#niter.goal = 5, 
				#
				stopRule = FALSE, 
				stopValue = 0.05, 
				maxAttempts = Inf
				#
				)
			}
		doRun_out[[i]] <- doRun_out_ThisTree	
		}
	#
	###########################################################
	#
	# record nDoRun as an attribute
	attr(doRun_out, "nDoRun") <- nDoRun
	#
	return(doRun_out)
	}
	