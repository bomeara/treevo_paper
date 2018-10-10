###################################################
# Individual Empirical Analyses and Simulations
##################################################

setwd("d://dave//workspace//treevo_paper//")

library(ape)
library(TreEvo)

######################################
# get empirical data
	
# 1) Anolis
	# repulsion - adaptive landscape dynamics - multi optima	

# obtain anolis tree - from Poe et al. 2017 (SystBiol)
	# their time-tree
anolisTree<-read.tree(file="datasets//anolis_PoeEtAl2018_datedMCCw0.5burnin.tre")

# obtain anolis trait data - 
	# Snout-Vent body-size data from Poe et al. 2018 (AmNat)
anolisTrait<-read.table("datasets//anolis_lntraits_matched_tabdelim_07-24-18.txt",
	header=TRUE,row.names=1)
anolisSize<-anolisTrait[,1]

# 2) Aquilegia
	# whittall et al. model of nectar spur increase in size

# obtain aquilegia tree (from Whittall and Hodges 2007?)
aquilegiaTree<-read.tree("datasets//aquilegia_Whttall&Hodges2007_figuredMCC.tre")

# obtain aquilegia trait data (from Whittall and Hodges 2007?)
	# need both nectur spur lengths and regime data
	#
aquilegiaTrait<-read.table("aquilegia_traitData.txt", header=FALSE, row.names=1)

# get just nectur spur length
aquilegiaSpurLength<-aquilegiaTrait[,2]
# and take the natural log
	# (note that the third column of the table was already the natural log)
	# previous code from Brian had 'log(data[,3])' - log of a log
aquilegiaSpurLength<-log(aquilegiaSpurLength)

# aquilegia regimes - pollinator syndromes
aquilegiaPollinators<-aquilegiaTrait[,14]
# regimes coded 0, 1, 2
	# 0 is bumble-bee, 1 is humming-bird, 2 is hawkmoth
# this probably won't be used?
# could use for post-analysis comparisons? Hmm
	
	
######	
# old aquilegia code
#assume generation time of 10 years (its a perennial plant), 
	# following Cooper et al. Plos ONe 2010 
	# Genetic Variation at Nuclear loci fails to distinguish group is about 3 MY,
		# phy height is 3. So each unit = 1,000,000 years or thus 100,000 generations
# TreeYears=100000
# timeStep<-1/TreeYears
# totalTreeLength=TreeYears*sum(phy$edge.length) #how many generations are represented

# number of expected polinator shifts based on parsimony is 7:
# parsimonyShifts=7
# pollinatorShiftRate=parsimonyShifts/totalTreeLength

	

	

#####	
	
##############################################################################

# need to reconstruct regimes down the aquilegia tree







###############################################################################
# generate sets of ideal trees for doing simulations on
   # idealTreeSets = c("Ideal-Balanced", "Ideal-Pectinate", "Ideal-Star") 
   # nTipSets = c(8, 16, 64)   
   
idealTrees<-list(
    #
    balanced_n8 = stree(n=8, type = "balanced", tip.label = NULL),
    balanced_n16 = stree(n=16, type = "balanced", tip.label = NULL),
    balanced_n64 = stree(n=64, type = "balanced", tip.label = NULL),
    #
    pectinate_n8 = stree(n=8, type = "left", tip.label = NULL),
    pectinate_n16 = stree(n=16, type = "left", tip.label = NULL),
    pectinate_n64 = stree(n=64, type = "left", tip.label = NULL), 
    #
    star_n8 = stree(n=8, type = "star", tip.label = NULL),
    star_n16 = stree(n=16, type = "star", tip.label = NULL),
    star_n64 = stree(n=64, type = "star", tip.label = NULL)
    )
	
	
	
	
	
	
	
	
	
	
### tree sets
# unresolved question? tree sizes
	# 10, 25, 50, 100, 200?, 500?
	# for now let's do first four
	# need to be power of two for balanced trees with stree
idealTreeSizes<-c(10,20,50,100)
nIdealTreeSet<-20
# should these all be extant only trees? Hmm.
makeExtantSTree<-function(tree){
	tree@edge.length<-runif(Nedge(tree))
	
	}


testTree<-list()
for(i in 1:length(idealTreeSizes)){
	treesOneSize<-list()
	# generate idealized trees
	# idealized - balanced tree
	treesOneSize[[1]]<-lapply(1:nIdealTreeSet,function(x) 
		stree(idealTreeSizes,type="balanced"))
	
		
	# idealized - pectinate tree
	treesOneSize[[1]]<-lapply(1:nIdealTreeSet,function(x) 
		stree(idealTreeSizes,type="left"))
	
	# idealized - star tree
	
	#
	# convert to multiPhylo, compress
	for(j in 1:length(treesOneSize)){
		class(treesOneSize[[j]])<-"multiPhylo"
		treesOneSize[[j]]<-.compressTipLabel(treesOneSize[[j]])
		
		}
	# 
	}
		
	
	
	
	
	
	
	
	
# realistic tree set
	# empirical neo tree
	# empirical paleo tree?

	# multiple empirical trees of varying sizes?




	
	
	
###################################################
###### TWO tests of treevo
#################################################

# Unresolved question: Number of particles? Number of generations?

#########################
### (1) test basic BM
	# fixed sigma square (rate)

	#Two different BM priors: unif with truth at say 25th percentile, rexp with mean not at true value (just to be realistic). Similar for root state.

	
#############################	
### (2) interaction model - repulsion between species with max bound 
	# assume log transformed traits, so no min bound
	# use realistic tree set

# Some models with 
	# no actual repulsion
	# moderate (species cross in trait space) 
	# high (no touching happens). 
# Models vary distance of max to the trait data
	# max is very close (most species bounce off it based on starting values), 
	# max is moderately far (start hitting in last 25% of sim), 
	# very far (never hit)
# That’s nine different parameter values, 
	# all of which I guess are scaled effective to the BM rate. 
	# Could do it as six: 
		#1) moderate repulsion, max bound is very close
		#2) moderate repulsion, max bound is moderately close
		#3) moderate repulsion, max bound is very far
		#4) no repulsion, max bound is moderately close
		#5) high repulsion, max bound is moderately close 		
		
#########################
# 3) Time - autoregressive model with optimum based on a factor that changes through time (like O2 concentration)
 # empirical tree, like in simulations

 # Parameters of the regression: 
	# a single variable function to convert O2 to optimal gill size or whatever
	# strength of pull
	# BM wiggle
	
# abd presumably the tracked env factor will be analyzed many times over





	
	
	

#' @rdname intrinsicModels
#' @export
autoregressiveMultOPtimaIntrinsic <- function(params, states, timefrompresent) {
    #a discrete time OU, same sd, mean, and attraction for all chars
    #params[1] is sd (sigma), params[2] is attractor (ie. character mean), params[3] is attraction (ie. alpha)
    sd <- params[1]
    attractor <- params[2]
    attraction <- params[3]    #in this model, this should be between zero and one
    newdisplacement <- rpgm::rpgm.rnorm(n = length(states), mean = (attractor-states)*attraction, sd = sd) #subtract current states because we want displacement
    return(newdisplacement)
    }   
   	
	
	
	
