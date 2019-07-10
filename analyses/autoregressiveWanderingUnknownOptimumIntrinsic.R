
R version 3.6.0 (2019-04-26) -- "Planting of a Tree"
Copyright (C) 2019 The R Foundation for Statistical Computing
Platform: x86_64-pc-linux-gnu (64-bit)

R is free software and comes with ABSOLUTELY NO WARRANTY.
You are welcome to redistribute it under certain conditions.
Type 'license()' or 'licence()' for distribution details.

  Natural language support but running in an English locale

R is a collaborative project with many contributors.
Type 'contributors()' for more information and
'citation()' on how to cite R or R packages in publications.

Type 'demo()' for some demos, 'help()' for on-line help, or
'help.start()' for an HTML browser interface to help.
Type 'q()' to quit R.

> ###################################################
> # Individual Empirical Analyses and Simulations
> ##################################################
> 
> # Control Box
> 
> # number of simulated trait datasets to do for simulated-trait runs
> nSimTrait <- 10
> 
> # error for run with mis-specified prior on sigmasq in a pure-BM model
> 	# the mean of the normal prior is multiplied by this value
> 	# 100 = mean of rate prior is off by two orders of magnitude!
> ratePriorError <- 100
> 
> # root age for idealized simulated trees from time=0
> 	# similar to Anolis root depth (51.49056)
> idealTreeDepth <- 50    
> 
> # simulation resolution
> 	# recc default is 1000
> generation.time <- 200000 
> 
> # control parameters for multicore and simulation resolution
> multicore <- TRUE 
> coreLimit <- 6
> 
> # control parameters for MCMC / ABC
> nRuns <- 1            		  # use 2 - recc default is 2 
> 								#(for testing, use 1)
> nStepsPRC <- 2        		  # use 5 - recc default is 5 
> 								#(for testing, use 2)
> numParticles <- 5 			  # use 300 - recc default is 300 
> 								#(for testing, use 5)
> nInitialSimsPerParam <- 10 	  # use 100 - recc default is 100 
> 								#(for testing, use 10)
> nInitialSims <- 5			  # use NULL - default is NULL = 100 per param 
> 								#(for testing, use 5)
> 
> #### miscellaneous controls
> # save data during runs?
> saveData <- FALSE
> # print out progress to terminal?
> verboseParticles <- FALSE
> 
> ###########################
> # FOR CONTINUING FROM A PREVIOUS TEST
> # (...if output from previous analyses exist at all)
> continueFromPrevious <- FALSE
> 
> source(
+ 	"~/treevo_paper/analyses/simulations_framework_script.R"
+ 	)
Error in file(filename, "r", encoding = encoding) : 
  cannot open the connection
Calls: source -> file
In addition: Warning message:
In file(filename, "r", encoding = encoding) :
  cannot open file '/home/bomeara/treevo_paper/analyses/simulations_framework_script.R': No such file or directory
Execution halted
