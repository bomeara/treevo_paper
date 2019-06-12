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

# simulation resolution
	# recc default is 1000
generation.time <- 200000 

# control parameters for multicore and simulation resolution
multicore <- TRUE 
coreLimit <- 6

# control parameters for MCMC / ABC
nRuns <- 1            		# recc default is 2
nStepsPRC <- 2        		# recc default is 5
numParticles <- 5 			# recc default is 300
nInitialSimsPerParam <- 10 	# recc default is 100
nInitialSims <- 5			# recc default is NULL (100 per param)

#### miscellaneous controls
# save data during runs?
saveData <- FALSE
# print out progress to terminal?
verboseParticles <- FALSE

source(
	"d://dave//workspace//treevo_paper//analyses//simulations_framework_script.R"
	)
