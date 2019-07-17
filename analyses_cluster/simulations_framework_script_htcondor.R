# simulation framework
# (shouldn't need to be changed)

##################################################
# Control Box

# number of simulated trait datasets to do for simulated-trait runs
nSimTrait <- 10

# error for run with mis-specified prior on sigmasq in a pure-BM model
# the mean of the normal prior is multiplied by this value
# 100 = mean of rate prior is off by two orders of magnitude!
ratePriorError <- 100

# root age for idealized simulated trees from time=0
# similar to Anolis root depth (51.49056)
idealTreeDepth <- 50

# simulation resolution
# recc default is 1000
generation.time <- 1000

# control parameters for multicore and simulation resolution
multicore <- TRUE
coreLimit <- 24

# control parameters for MCMC / ABC
nRuns <- 2            		  # use 2 - recc default is 2
#(for testing, use 1)
nStepsPRC <- 5        		  # use 5 - recc default is 5
#(for testing, use 2)
numParticles <- 300 			  # use 300 - recc default is 300
#(for testing, use 5)
nInitialSimsPerParam <- 100 	  # use 100 - recc default is 100
#(for testing, use 10)
nInitialSims <- 100			  # use NULL - default is NULL = 100 per param
#(for testing, use 5)

#### miscellaneous controls
# save data during runs?
saveData <- FALSE
# print out progress to terminal?
verboseParticles <- FALSE

###########################
# FOR CONTINUING FROM A PREVIOUS TEST
# (...if output from previous analyses exist at all)
continueFromPrevious <- FALSE


library(ape)
library(TreEvo)

# get package versions
if(packageVersion("TreEvo") < "0.21.0"){
	stop("Update TreEvo first!")
}

message(paste0(
	"TreEvo Version Used: ",
	packageVersion("TreEvo")
))
message(paste0(
	"ape Version Used: ",
	packageVersion("ape")
))

#setwd("/share/bomeara/treevo_paper//")
source("functions_for_analysis.R")
source("functions_for_aquilegia_models.R")

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
	file="..//datasets//anolis_PoeEtAl2018_datedMCCw0.5burnin.tre"
)
#
# obtain anolis trait data -
# Snout-Vent body-size data from Poe et al. 2018 (AmNat)
anolisTrait <- read.table(
	"../datasets//anolis_lntraits_matched_tabdelim_07-24-18.txt",
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
	message(paste0(
		"The following OTUs (",
		length(droppers),
		") on the Anolis tree do not appear to\n",
		" have size data and thus will be dropped: \n ",
		paste0(strwrap(
			paste0(droppers, collapse=", ")
		),collapse="\n ")
	))
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
	"../datasets//aquilegia_Whttall&Hodges2007_figuredMCC.tre"
)
# need to clear away the root edge length
aquilegiaTree$root.edge <- NULL
#
# obtain aquilegia trait data (from Whittall and Hodges 2007?)
# need both nectur spur lengths and regime data
#
aquilegiaTrait <- read.table(
	"../datasets//aquilegia_traitData.txt",
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
	message(paste0(
		"The following OTUs(",
		length(droppers),
		") on the Aquilegia tree do not appear to\n",
		" have spur length data and thus will be dropped: \n ",
		paste0(strwrap(
			paste0(droppers, collapse=", ")
		),collapse="\n ")
	))
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
#
idealTrees <- list(
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
	),
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
	)
)
#
# all of these need to have edge lengths
idealTrees <- lapply(idealTrees, compute.brlen)
# multiple edge lengths by 50
idealTrees <- lapply(idealTrees,
	function(x) {
		x$edge.length <- x$edge.length * idealTreeDepth
		return(x)
	}
)
#
# make all trees artificially bifurcating
idealTrees <- lapply(idealTrees, multi2di)
#
# test that they are ultrametric
if(!all(sapply(idealTrees,is.ultrametric))){
	stop("Not all idealized simulated trees came out as ultrametric ?!")
}
#
# make multiPhylo
class(idealTrees) <- "multiPhylo"
#
# compress tip labels? No, I don't think that works for trees of different sizes
#	trees <- .compressTipLabel(trees)
#
######################################################################################
message("##############################")
message("######### Beginning Analyses ############")
#
# time to get table, process the inputs listed
#
# get simulation run table
simRunTable <- read.csv(
	file="simulation_sets_parameters_table.csv",
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
# start an empty analysisSetup
analysisSetup <- list()
loadAnalysisSetup <- FALSE
#
#
############################################
#
#
# if(continueFromPrevious){
# 	#
# 	# output
# 	outFiles <- file.info(list.files(
# 		"./saved_output/", full.names = TRUE
# 		))
# 	outFiles <- rownames(outFiles)[which.max(outFiles$mtime)]
# 	# if there are any files...
# 	if(length(outFiles) > 0){
# 		# replace analysisOutput
# 		analysisOutputOld <- readRDS(file=outFiles)
# 		if(analysisOutputOld == "analysisOutput"){
# 			stop("Somehow the old output is just the string 'analysisOutput'...?")
# 			}
# 		if((length(analysisOutputOld) == length(analysesNames))
# 				& identical(analysesNames,names(analysisOutput))){
# 			message("Loading output file from previous run...")
# 			analysisOutput <- analysisOutputOld
# 			#
# 			# also load old analysisSetup
# 			loadAnalysisSetup <- TRUE
# 		}else{
# 			warning(paste0(
# 				"Format of previous output file does not match current script expectations\n",
# 				"Beginning analyses without loading previous output file..."
# 				))
# 			if(any(sapply(analysisOutput,length) > 1)){
# 				stop("How did analysisOutput get non-fresh data without restarting from old?")
# 				}
# 			}
# 	}else{
# 		message(paste0(
# 			"No output files from previous runs found\n",
# 			"Beginning analyses without loading any previous output file..."
# 			))
# 		}
# 	}
#
#
# also load old analysisSetup if loading old output
# if(loadAnalysisSetup){
# 	# output
# 	outFiles <- file.info(list.files(
# 		"./saved_setup/", full.names = TRUE
# 		))
# 	outFiles <- rownames(outFiles)[which.max(outFiles$mtime)]
# 	# if there are any files...
# 	if(length(outFiles) > 0){
# 		# replace analysisOutput
# 		analysisSetup <- readRDS(file=outFiles)
# 		#
# 		message("Loading analysis setup from previous run...")
# 	}else{
# 		message("No previous analysis setup found.")
# 		}
# 	}
# #
# ################################################################
# # test that analysis output is useable
# if(!identical(analysesNames, names(analysisOutput))){
# 	stop(
# 		"analysisOutput seems to be corrupt - names do not match analysesNames"
# 		)
# 	}
# #
##############################################
# make a new save file name for output
# saveFileName <- paste0(
# 	".//saved_output//",
# 	"analysisOutput_saved_",
# 	format(Sys.time(), "%m-%d-%y"),
# 	".rds"
# 	)
# #
# # save initial file
# saveRDS(analysisOutput,
# 	file = saveFileName
# 	)
# #
# # do the same for the setup file
# saveSetupName <- paste0(
# 	".//saved_setup//",
# 	"analysisSetup_saved_",
# 	format(Sys.time(), "%m-%d-%y"),
# 	".rds"
# 	)
# #
# # save initial file
# saveRDS(analysisSetup,
# 	file = saveSetupName
# 	)
#
######################################
#
# Let's runs the analyses!
#
# run all independent analyses
message("###############################################")
message("#########  Independent Analyses  ##############")
#
analysisSetup <- c()
for (i in whichIndependentPrevRun){
	if(analysisOutput[[i]] == analysesNames[i]){
		#
		message("#######################################")
		message("######   Now running -- ", analysesNames[i], "  #########")
		#
		runParameters <- simRunTable[i, , drop = FALSE]
		#
			analysisSetup <- setupRunAnalysis(
				runParameters = runParameters,
				#
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
				indepAnalyses_extrinsicOut = NULL
			)
			#
			# save analysisSetup
			# saveRDS(analysisSetup,
			# 	file = saveSetupName
			# )

		#################
		#
		# description of framework for sending jobs to htcondor
		# taken from gitter chat with Brian
		#
		#
		# So, line 422 of the htcondor.R file
		# I'm doing your loop, but always doing a new setup file 
			# I save this at line 456
		# line 458, I write a new R script 
			# which ends at line 500
		# line 502, I make a new bash script
		# line 507, I make it executable
		# line 502-520, I make a condor submit script
		# and line 522, submit the job
			# so, the submit script has to call a program. 
			# That's the sh script, which just says to run Rscript with the given arguments
			# The submit script has to know what files to copy over and those files must exist
			# thus the line 500 file and the line 456 file
		# so condor submits the submit script, which says, "hey, I need 24 cores" 
			# and gets assigned a machine, then it copies over the files, 
			# calls the bash script, which calls Rscript, which then runs. 
			# When done (or dead), condor takes any files the script has made and 
			# copies them over to the original directory		

		# I'm doing your loop, but always doing a new setup file 
			# I save this at line 456
		save(list=ls(), 
			file=paste0(
				"Data_",
				analysesNames[i],
				"_",
				Sys.Date(),
				".rda"
				)
			)

		###################################################
		# make the R script
		#
		# line 458, I write a new R script 
			# which ends at line 500
		#
		cat(
			paste0(
#############################################		
'
library(ape)
library(TreEvo)
load("Data_',
	analysesNames[i],
	"_",
	Sys.Date(),
	'.rda")

# get package versions
if(packageVersion("TreEvo") < "0.21.0"){
	stop("Update TreEvo first!")
}

message(paste0(
	"TreEvo Version Used: ",
	packageVersion("TreEvo")
))
message(paste0(
	"ape Version Used: ",
	packageVersion("ape")
))

# now doRun!
result <- doRunAnalysis(
	treeList = analysisSetup$treeList,
	traitDataList = analysisSetup$traitDataList,
	runLabel = analysisSetup$runLabel,
	nDoRun = analysisSetup$nDoRun,
	intrinsicFunctionToFit = analysisSetup$intrinsicFunctionToFit,
	extrinsicFunctionToFit = analysisSetup$extrinsicFunctionToFit,
	intrinsicArgList = analysisSetup$intrinsicArgList,
	extrinsicArgList = analysisSetup$extrinsicArgList,
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
save(result, file="Results_', 
	analysesNames[i],
	"_",
	Sys.Date(), 
	'.rda")
'
#######################################
				), 
			file = paste0(
				"Run_",
				analysesNames[i],
				"_",
				Sys.Date(),
				".R"
			))

		###########################################
		# make the sh file 
			# line 502, I make a new bash script
			#
			# That's the sh script, which just says to run Rscript with the given arguments
			#
		#########	
		cat(
			paste0(
##################################		
'#!/bin/bash

Rscript Run_',
	analysesNames[i],
	"_",
	Sys.Date(),
	'.R'
#####################################
				), 
			file=paste0(
				"Run_",
				analysesNames[i],
				"_",
				Sys.Date(),
				".sh"
				))
				
		##############################################
		#
		
		
		
		# line 502-520, I make a condor submit script		

# line 507, I make the sh file executable

		system(paste0(
			"chmod u+x Run_",
			analysesNames[i],
			"_",
			Sys.Date(),
			".sh"
			))
			
		#####################################################
		# Make a qsub file - SUBMIT SCRIPT
		#
		# The submit script has to know what files to copy over and those files must exist
			# thus specifying the line 500 file and the line 456 file		
		#
		cat(
			paste0(
#####################################################
'log = run_$(Cluster)_$(Process)_', 
	analysesNames[i],
	"_",
	Sys.Date(), 
	'.log/n',
'error = run_$(Cluster)_$(Process)_', 
	analysesNames[i],
	"_",
	Sys.Date(), 
	'.err/n',
'output = run_$(Cluster)_$(Process)_', 
	analysesNames[i],
	"_",
	Sys.Date(), 
	'.out/n',
'\n,'
'executable = Run_',
	analysesNames[i],
	"_",
	Sys.Date(),
	'.sh/n',
'
should_transfer_files = YES
when_to_transfer_output = ON_EXIT/n',
'
transfer_input_files = Data_',
	analysesNames[i],
	"_",
	Sys.Date(),
	'.rda,Run_',
	analysesNames[i],
	"_",
	Sys.Date(),
	'.R,Run_',
	analysesNames[i],
	"_",
	Sys.Date(),
	'.sh/n',
'
request_cpus = 24
queue 1
'
###############################################
				), 
			file=paste0(
				"Run_",
				analysesNames[i],
				"_",
				Sys.Date(),
				".qsub"
				))

		#####################################
		# SUBMITTING THE JOB
		#
		# and line 522, submit the job
			# so, the submit script has to call a program. 
			# That's the sh script, which just says to run Rscript with the given arguments
			# The submit script has to know what files to copy over and those files must exist
			# thus the line 500 file and the line 456 file		
		# so condor submits the submit script, which says, "hey, I need 24 cores" 
			# and gets assigned a machine, then it copies over the files, 
			# calls the bash script, which calls Rscript, which then runs. 
			# When done (or dead), condor takes any files the script has made and 
			# copies them over to the original directory		
		######
		system(paste0(
			"/usr/bin/condor_submit Run_",
			analysesNames[i],
			"_",
			Sys.Date(),
			".qsub"
			))
		}
	}
