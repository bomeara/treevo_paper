library(ape)
		library(TreEvo)
		load("Data_An_Emp_Bound_2019-07-14.rda")

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
		save(result, file="Results_An_Emp_Bound_2019-07-14.rda")
		