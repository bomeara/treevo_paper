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


