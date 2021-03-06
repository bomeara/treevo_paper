# This function calculates Effective Sample Size (ESS) on results.  
# Performs the best when results are from multiple runs.

pairwiseESS(results$particleDataFrame)

# bayesCoverageProb
	# for comparing true / generating parameter values to the posteriors of analyses done on that data

# plotUnivariatePosteriorVsPrior 
	# plot priors versus their posteriors - useful for runs with bad prior on BM?
	
```
# examples of getting density coordinates and summary statistics from distributions
priorKernal<-getUnivariatePriorCurve(priorFn="normal", priorVariables=c(28,2),
	nPoints=100000, from=NULL, to=NULL, prob=0.95)
postKernal<-getUnivariatePosteriorCurve(acceptedValues=results$particleDataFrame$starting_1,
	from=NULL, to=NULL, prob=0.95)
# let's compare this (supposed) prior against the posterior in a plot
plotUnivariatePosteriorVsPrior(posteriorCurve=postKernal, priorCurve=priorKernal,
	label="parameter", trueValue=NULL, prob=0.95)
```
	
# plotPosteriors
	# for each free parameter in the posterior, a plot is made of the distribution of values estimate in the last generation
	# can also be used to visually compare against true (generating) parameter values in a simulation.
	
plotPosteriors(particleDataFrame=resultsBM$particleDataFrame,
   priorsMat=resultsBM$PriorMatrix)

# highestPostDens 
	# get weighted mean, standard deviation, upper and lower highest posterior density (HPD) for each free parameter in posterior. 

highestPostDens(results$particleDataFrame, percent=0.95, returnData=FALSE)
	
# plotABC_3D
	# Plot posterior density distribution for each generation in 3d plot window 


#############################################################################################
# notes from conversation with Peter Smits (05-09-18)
#
# so I'm doing approximate bayesian computation and the question is, what do I want to show to the reader
# I want to show posterior parameter estimates from real data, 
	# and show that they are very different from parameter estimates made under other models,
	# or under the same model but with simulated data, for scenarios with a small number of models
# what i do is make the same series of posterior predictive checks
	# and demonstrate how your prefered model better recapitulates the data it was fit to
#
# 
############################################################
# ECDF
# ECDF - empirical cumulative distribution function = the ranked order accumulation curve
# http://stat.ethz.ch/R-manual/R-devel/library/stats/html/ecdf.html
# ecdf is a cool way of summarizing the entire dataset graphically
#
# how well does simulations under a fit model reproduce ecdf or the density of the original data? 
# if your model does a better job of doing that, then it is straight up a better model
# it also goes beyond parameter estimates and towards the model describing the data
#
# bayes wants to describe more than just the expected value. it is greedy and wants to describe the whole posterior
# the posterior predictive distribution describes all data sets that are consistent with the model, given the original input information
# if the PPD doesn't look like the empirical data, then the model is not describing your data
#####################################################################################
# from 06-21-18
# okay so the general sketch is particles from the posterior, simulate under this set of parameters N times, 
    # and compare the original ECDF for each parameter to the simulated
# my other idea:
# draw parameter estimates from some posterior particle, simulate under those parameters, 
    #  then test if 'true' generating parameters are actually within the 95% HDF of the simulated posterior
    # deals with how we don't really understand how adequate the models are for giving unbiased estimates of parameters

# checkAdequacy              # sixAnalysesTwoModels
    

#checkAdequacy <- function(){
#    }



#I mean, writing a function that just takes arguments: tree, params, etc. and returns results would be good

#a lot of this is (dataset) and six corresponding analyses

#fit model A, model B to real data, then simulate under model A, model B and fits both model A and B to both
#(where A is usually BM)

