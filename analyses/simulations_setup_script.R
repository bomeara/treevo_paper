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

# root age for idealized simulated trees from time=0
	# similar to Anolis root depth (51.49056)
idealTreeDepth <- 50    

# simulation resolution
	# recc default is 1000
generation.time <- 200000 

# control parameters for multicore and simulation resolution
multicore <- TRUE 
coreLimit <- 6

# control parameters for MCMC / ABC
nRuns <- 1            		  # use 2 - recc default is 2 (for testing, use 1)
nStepsPRC <- 2        		  # use 5 - recc default is 5 (for testing, use 2)
numParticles <- 5 			  # use 300 - recc default is 300 (for testing, use 5)
nInitialSimsPerParam <- 10 	  # use 100 - recc default is 100 (for testing, use 10)
nInitialSims <- 5			  # use NULL - recc default is NULL-ie 100 per param (for testing, use 5)

#### miscellaneous controls
# save data during runs?
saveData <- FALSE
# print out progress to terminal?
verboseParticles <- FALSE

###########################
# FOR CONTINUING FROM A PREVIOUS TEST
# (...if output from previous analyses exist at all)
continueFromPrevious <- TRUE

source(
	"~/treevo_paper/analyses/simulations_framework_script.R"
	)
