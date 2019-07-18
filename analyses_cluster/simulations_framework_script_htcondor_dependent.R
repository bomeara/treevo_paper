

##################################################
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

setwd("/share/bomeara/treevo_paper//")

############################################

# load everything 




files <- list.files()


# find all the old workspaces

# find the original saved workspace

# load things to get analysisOutput
	# a list that contains all dependent run data







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
#
# BUT NOTICE THAT SOME OF THESE DO NOT HAVE 
	# CORRESPONDING INDEP ANALYSES !
#
# actual indep analyses performed:
	# An_Emp_DispBound
	# An_Emp_Bound
	# An_Emp_BrownMotion
	# An_Emp_Disp
	# An_Emp_TimeReg
	# Aq_Emp_3Opt2Bound
	# Aq_Emp_BrownMotion
#
# ones not covered by indep analyses
	# An_Emp_Bound_BoundByStartingState
	# An_Emp_Bound_BoundByMinValue
	# An_Emp_Bound_BoundOneRangeAway
#
#############################
# get the stuff necessary for doing the dependent analyses
#
# get model parameters from runs 
	# that will be used for dependent simulations
# use extract on all indep analyses now
	# then can call these later for dependent analyses
	# without have to extract same data many times
#
# Note that following functions will only look at first analysis
	# this doesn't matter - all indep analyses should only have one analysis
	# one empirical tree, one empirical trait, thus only one analysis to examine
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
#
# make sure named correctly
names(indepAnalyses_intrinsicOut) <- analysesNames[whichIndependentPrevRun]
names(indepAnalyses_extrinsicOut) <- analysesNames[whichIndependentPrevRun]
#
# add intrinsic models not included
	# An_Emp_Bound_BoundByStartingState
	# An_Emp_Bound_BoundByMinValue
	# An_Emp_Bound_BoundOneRangeAway
#
# all of these are based on An_Emp_Bound
boundInt <- indepAnalyses_intrinsicOut$An_Emp_Bound
#
# An_Emp_Bound_BoundByStartingState
	# bound is right by the starting state, leading to
	# diffusion away from left-hand wall dynamics
boundIntStarting <- boundInt	
# set bound equal to starting state
boundIntStarting$intrinsicValues['intrinsic_2'] <- boundIntStarting$startingValues[1]
indepAnalyses_intrinsicOut$An_Emp_Bound_BoundByStartingState <- boundIntStarting
#
# An_Emp_Bound_BoundByMinValue
	# bound is at the minimum value observed for anolisSize
boundIntMin <- boundInt	
# set bound equal to minimum size observed
boundIntMin$intrinsicValues['intrinsic_2'] <- min(anolisSize)
indepAnalyses_intrinsicOut$An_Emp_Bound_BoundByMinValue <- boundIntMin
#
# An_Emp_Bound_BoundOneRangeAway
	# what if the bound was very distant -
		# i.e. one range (max-min) away from the min
oneRange <- max(anolisSize) - min(anolisSize)
oneRangeAway <- min(anolisSize) - oneRange
boundOneR <- boundInt	
# set bound equal to minimum size observed
boundOneR$intrinsicValues['intrinsic_2'] <- oneRangeAway
indepAnalyses_intrinsicOut$An_Emp_Bound_BoundOneRangeAway<- boundOneR
#
################################################################

# run all dependent analyses
#
message("#############################################")
message("#########  Dependent Analyses  ##############")
#
for (i in whichDependentPrevRun){
	if(analysisOutput[[i]] == analysesNames[i]){
		#
		message("#####################################")
		message("######   Now running -- ", analysesNames[i], "  ##########")
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
			indepAnalyses_intrinsicOut = 
				indepAnalyses_intrinsicOut,
			indepAnalyses_extrinsicOut = 
				indepAnalyses_extrinsicOut
			)	
			


		#################
		# now doRun! 
		#
		analysisOutput[[i]] <-	doRunAnalysis(
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
		



		}
	}




			##################################################################################
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
			#
			######################################################################
			#
			# get file names
			# 	
			analysisName <- analysesNames[i]
			#
			# core file name
			run_file_name_generic <- paste0(
				analysesNames[i],
				"_",
				Sys.Date()
				)
			#
			# RDA file
			run_workspace_file <- paste0(
				"Data_",
				run_file_name_generic,
				".rda"
				)
			# file for saved results
			run_saved_results_file <- paste0(
				"Results_", 
				run_file_name_generic, 
				".rda"
				)
			# R file
			run_R_script <- paste0(
				"Run_",
				run_file_name_generic,
				".R"
				)
			# sh file
			run_sh_script <- paste0(
				'Run_',
				run_file_name_generic,
				'.sh'
				)	
			# condor submission script (qsub)
			run_submission_script <- paste0(
					"Run_",
					run_file_name_generic,
					".qsub"
					)
			# generic condor file name
			run_condor_generic_name <- paste0(
				'run_$(Cluster)_$(Process)_',
				run_file_name_generic
				)
			# error file
			run_error_condor <- paste0(
				run_condor_generic_name, 
				'.err'
				)
			# log file
			run_log_condor <- paste0(
				run_condor_generic_name, 
				'.log'
				)
			# output file
			run_output_condor <- paste0(
				run_condor_generic_name, 
				'.out'
				)
		#####################################################################
		# I'm doing your loop, but always doing a new setup file 
			# I save this at line 456
		save(list=ls(), file= run_workspace_file)
		#
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
load(run_workspace_file_name)

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
save(result, file=',run_saved_results_file,')'
####################################################
				), 
			file = run_R_script
			)

		###########################################################
		# make the sh file 
			# line 502, I make a new bash script
			#
			# That's the sh script,
				# which just says to run Rscript with the given arguments
			#
		#########	
		cat(
			paste0(
				"#!/bin/bash\n,"
				"Rscript ",run_R_script
				), 
			file = run_sh_script
			)
		#		
		##############################################
		#
		# line 507, I make the sh file executable
			# you mean you execute the sh file?
		#
		system(paste0("chmod u+x",run_sh_script))
		#	
		#####################################################
		# Make a qsub file - SUBMIT SCRIPT
		#
		# line 502-520, I make a condor submit script		
		#
		# The submit script has to know what files to copy over and those files must exist
			# thus specifying the line 500 file and the line 456 file		
		#
		cat(
			paste0(
#####################################################
'log = ',run_log_condor,'
error = ',run_error_condor,'
output = ', run_output_condor,'
executable = ', run_sh_script,'
should_transfer_files = YES
when_to_transfer_output = ON_EXIT,
transfer_input_files = ',
	run_workspace_file,",",
	run_R_script,",",
	run_sh_script,'
request_cpus = 24
queue 1
'
###############################################
				), 
			file=run_submission_script
			)

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
			"/usr/bin/condor_submit ",
			run_submission_script
			))
		}
	}


