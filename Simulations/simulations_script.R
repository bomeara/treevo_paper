# simulation and analyses 04-26-18

# load packages
library(ape)
library(TreEvo)

### tree sets
# unresolved question? tree sizes
	# 10, 25, 50, 100, 200?, 500?
	# for now let's do first four
	# need to be power of two for balanced trees with stree
idealTreeSizes<-c(10,20,50,100)
nIdealTreeSet<-20
# should these all be extant only trees? Hmm.
makeExtant<-


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

#####################################################
#############################################################
####################################################
# two empirical examples
	
# 1) Anolis
	# repulsion
	
# 2) Aquilegia
	# whittall et al. model of nectar spur increase in size


