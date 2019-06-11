
##################################################

library(ape)
library(TreEvo)
setwd("d://dave//workspace//treevo_paper//")
source(".//analyses//functions_for_analysis.R")

######################################
# get empirical data
#############################################
#	
# 1) Anolis
	# repulsion - adaptive landscape dynamics - multi optima	
#
# obtain anolis tree - from Poe et al. 2017 (SystBiol)
	# their time-tree
anolisTree <- read.tree(
	file="datasets//anolis_PoeEtAl2018_datedMCCw0.5burnin.tre"
	)
#
# obtain anolis trait data - 
	# Snout-Vent body-size data from Poe et al. 2018 (AmNat)
anolisTrait <- read.table(
	"datasets//anolis_lntraits_matched_tabdelim_07-24-18.txt",
	header=TRUE,row.names=1
	)
#
anolisSize <- anolisTrait[,1]             # ,drop = FALSE]
names(anolisSize) <- rownames(anolisTrait)
#
# need to remove all unshared taxa from the tree
#
# crop traits down to those in the tree
anolisSize <- anolisSize[anolisTree$tip.label]
names(anolisSize) <- anolisTree$tip.label
# are any NA?
anyMatchesNA <- is.na(anolisSize)
if(any(anyMatchesNA)){	
	droppers <- names(anolisSize)[anyMatchesNA]
	message(
		"The following OTUs on the Anolis tree do not appear to\n",
		" have size data and thus will be dropped: \n",
		paste0(droppers, collapse=", "))
	anolisTree <- drop.tip(anolisTree, droppers)
	anolisSize <- anolisSize[anolisTree$tip.label]
	names(anolisSize) <- anolisTree$tip.label
	}
# make into a multiPhylo list
anolisTreeList <- list(anolisTree = anolisTree)
class(anolisTreeList) <- "multiPhylo"
#
################################################
#
# 2) Aquilegia
	# whittall et al. model of nectar spur increase in size
#
# obtain aquilegia tree (from Whittall and Hodges 2007?)
aquilegiaTree <- read.tree(
	"datasets//aquilegia_Whttall&Hodges2007_figuredMCC.tre"
	)
#
# obtain aquilegia trait data (from Whittall and Hodges 2007?)
	# need both nectur spur lengths and regime data
	#
aquilegiaTrait <- read.table(
	"datasets//aquilegia_traitData.txt", 
	header=FALSE, row.names=1
	)
#
# get just nectur spur length
aquilegiaSpurLength <- aquilegiaTrait[,2]           # , drop = FALSE]
names(aquilegiaSpurLength) <- rownames(aquilegiaTrait)
# and take the natural log
	# (note that the third column of the table was already the natural log)
	# previous code from Brian had 'log(data[,3])' - log of a log
aquilegiaSpurLength <- log(aquilegiaSpurLength)
# crop traits down to those in the tree
aquilegiaSpurLength <- aquilegiaSpurLength[aquilegiaTree$tip.label]
names(aquilegiaSpurLength) <- aquilegiaTree$tip.label
#
# aquilegia regimes - pollinator syndromes
aquilegiaPollinators <- aquilegiaTrait[,14]
names(aquilegiaPollinators) <- rownames(aquilegiaTrait)
# crop traits down to those in the tree
aquilegiaPollinators <- aquilegiaPollinators[aquilegiaTree$tip.label]
names(aquilegiaPollinators) <- aquilegiaTree$tip.label
#
# regimes coded 0, 1, 2
	# 0 is bumble-bee, 1 is humming-bird, 2 is hawkmoth
# this probably won't be used directly?
# could use for post-analysis comparisons? Hmm
#
# need to remove all unshared taxa from the tree
	# will do this ONLY relative to spur length vector!
#

# are any NA?
anyMatchesNA <- is.na(aquilegiaSpurLength)
if(any(anyMatchesNA)){	
	droppers <- names(aquilegiaSpurLength)[anyMatchesNA]
	message(
		"The following OTUs on the Aquilegia tree do not appear to\n",
		" have spur length data and thus will be dropped: \n",
		paste0(droppers, collapse=", "))
	aquilegiaTree <- drop.tip(aquilegiaTree, droppers)
	# and drop from trait data
	aquilegiaSpurLength <- aquilegiaSpurLength[aquilegiaTree$tip.label]
	aquilegiaPollinators <- aquilegiaPollinators[aquilegiaTree$tip.label]
	names(aquilegiaSpurLength) <- aquilegiaTree$tip.label
	names(aquilegiaPollinators) <- aquilegiaTree$tip.label
	}
# make into a multiPhylo list
aquilegiaTreeList <- list(aquilegiaTree = aquilegiaTree)
class(aquilegiaTreeList) <- "multiPhylo"
#
###############################################
# legacy aquilegia code from Brian O'Meara:
# 
# assume generation time of 10 years (its a perennial plant), 
	# following Cooper et al. Plos ONe 2010 
	# Genetic Variation at Nuclear loci fails to distinguish group is about 3 MY,
# So =>> phy height is 3. 
# Thus each unit = 1,000,000 years or 100,000 generations
#
# TreeYears=100000
# timeStep <- 1/TreeYears
# totalTreeLength=TreeYears*sum(phy$edge.length) #how many generations are represented
# number of expected polinator shifts based on parsimony is 7:
# parsimonyShifts=7
# pollinatorShiftRate=parsimonyShifts/totalTreeLength
#
#
###############################################################################
# generate sets of ideal trees for doing simulations on
   # idealTreeSets = c("Ideal-Balanced", "Ideal-Pectinate", "Ideal-Star") 
   # nTipSets = c(8, 16, 64)   
   
idealTrees <- list(
    #
    balanced_n8 = stree(
		n=8, 
		type = "balanced", tip.label = NULL
		),
    balanced_n16 = stree(
		n=16, 
		type = "balanced", tip.label = NULL
		),
    balanced_n64 = stree(
		n=64, 
		type = "balanced", tip.label = NULL
		),
    #
    pectinate_n8 = stree(
		n=8, 
		type = "left", tip.label = NULL
		),
    pectinate_n16 = stree(
		n=16, 
		type = "left", tip.label = NULL
		),
    pectinate_n64 = stree(
		n=64, 
		type = "left", tip.label = NULL
		), 
    #
    star_n8 = stree(
		n=8, 
		type = "star", tip.label = NULL
		),
    star_n16 = stree(
		n=16, 
		type = "star", tip.label = NULL
		),
    star_n64 = stree(
		n=64, 
		type = "star", tip.label = NULL
		)
    )
# make multiPhylo	
class(idealTrees) <- "multiPhylo"
#
# compress tip labels? No, I don't think that works for trees of different sizes
	#	trees <- .compressTipLabel(trees)
######################################################################################
# time to get table, process the inputs listed
# 		
# get simulation run table
simRunTable <- read.csv(
	file="analyses//simulation_sets_parameters_table.csv",
	header=TRUE,
	stringsAsFactors=FALSE
	)	
#
# number of analyses
nAnalyses <- nrow(simRunTable)
# 
# names of analyses
analysesNames <- simRunTable$runLabel	
#
# which analyses are independent or dependent
whichIndependentPrevRun <- which(
	!as.logical(simRunTable$dependentPrevRun)
	)	
whichDependentPrevRun <- which(
	as.logical(simRunTable$dependentPrevRun)
	)
#
# create list for saving analysis output 
analysisOutput <- as.list(analysesNames)
# why numbers? use actual names
	# names(analysisOutput) <- 1:nAnalyses
names(analysisOutput) <- analysesNames
#
# Let's runs the analyses! 
#
# run all independent analyses
for (i in whichIndependentPrevRun){
	analysisOutput[[i]] <- runAnalysis(
		runParameters = simRunTable[i, , drop = FALSE],
		# inputs needed from script	above
		nSimTrait = nSimTrait,
		ratePriorError = ratePriorError,
		#
		anolisTreeList = anolisTreeList,
		anolisSize = anolisSize,
		aquilegiaTreeList = aquilegiaTreeList,	
		aquilegiaSpurLength = aquilegiaSpurLength,
		idealTrees = idealTrees,
		#
		indepAnalyses_intrinsicOut = NULL,
		indepAnalyses_extrinsicOut = NULL,
		#
		# presets
		generation.time = generation.time,
		multicore = multicore,
		coreLimit = coreLimit,				
		nRuns = nRuns, 
		nStepsPRC = nStepsPRC, 
		numParticles = numParticles, 
		nInitialSimsPerParam = nInitialSimsPerParam, 
		nInitialSims = nInitialSims, 
		saveData = saveData,
		verboseParticles = verboseParticles

		)
	#
	save(analysisOutput, 
		file = paste0("//saved_output//",
			"analysisOutput_saved_",
			format(Sys.time(), "%m-%d-%y"),
			".Rdata")
		)
	}
#############################
# dependent analyses
########################################
#indep runs that dep runs depend on :
#
# INTRINSIC
	# An_Emp_BrownMotion
	# An_Emp_Disp
	# An_Emp_Bound
	# An_Emp_DispBound
	# An_Emp_Bound_BoundByStartingState
	# An_Emp_Bound_BoundByMinValue
	# An_Emp_Bound_BoundOneRangeAway
	# An_Emp_TimeReg
	# Aq_Emp_3Opt2Bound
	# Aq_Emp_BrownMotion
# EXTRINSIC
	# An_Emp_DispBound
	# An_Emp_Disp
#############################
# get the stuff necessary for doing the dependent analyses
#
# get model parameters from runs 
	# that will be used for dependent simulations
# use extract on all indep analyses now
	# then can call these later for dependent analyses
	# without have to extract same data many times
#
indepAnalyses_intrinsicOut <- lapply(
	analysisOutput[whichIndependentPrevRun],
	extractIntrinsic_from_prcOut
	)
#
indepAnalyses_extrinsicOut <- lapply(
	analysisOutput[whichIndependentPrevRun],
	extractExtrinsic_from_prcOut
	)
# make sure named correctly
names(indepAnalyses_intrinsicOut) <- analysesNames[whichIndependentPrevRun]
names(indepAnalyses_extrinsicOut) <- analysesNames[whichIndependentPrevRun]
#
# run all dependent analyses
for (i in whichDependentPrevRun){
	analysisOutput[[i]] <- runAnalysis(
		runParameters = simRunTable[i, , drop = FALSE],
		# inputs needed from script	above
		nSimTrait = nSimTrait,
		ratePriorError = ratePriorError,
		#
		anolisTreeList = anolisTreeList,
		anolisSize = anolisSize,
		aquilegiaTreeList = aquilegiaTreeList,	
		aquilegiaSpurLength = aquilegiaSpurLength,
		idealTrees = idealTrees,
		#
		indepAnalyses_intrinsicOut = indepAnalyses_intrinsicOut,
		indepAnalyses_extrinsicOut = indepAnalyses_extrinsicOut,
		#
		# presets
		generation.time = generation.time,
		multicore = multicore,
		coreLimit = coreLimit,				
		nRuns = nRuns, 
		nStepsPRC = nStepsPRC, 
		numParticles = numParticles, 
		nInitialSimsPerParam = nInitialSimsPerParam, 
		nInitialSims = nInitialSims,
		saveData = saveData,
		verboseParticles = verboseParticles
		)
	#
	save(analysisOutput, 
		file = paste0("//saved_output//",
			"analysisOutput_saved_",
			format(Sys.time(), "%m-%d-%y"),
			".Rdata")
		)
	}
