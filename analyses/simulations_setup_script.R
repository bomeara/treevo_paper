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
generation.time <- 100000 

# control parameters for multicore and simulation resolution
multicore <- TRUE 
coreLimit <- 6

# control parameters for MCMC / ABC
nRuns <- 2 
nStepsPRC <- 3 
numParticles <- 20 
nInitialSimsPerParam <- 10 
StartSims <- 10

source(
	"d://dave//workspace//treevo_paper//analyses//simulations_framework_script.R"
	)
