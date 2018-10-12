###################################################
# Individual Empirical Analyses and Simulations
##################################################

# Control Box

# number of simulated trait datasets to do for simulated-trait runs
nSimTrait <- 10


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
# make into a multiPhylo list
anolisTreeList <- list(anolisTree = anolisTree)
class(anolisTreeList) <- "multiPhylo"

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

# legacy aquilegia code from Brian O'Meara:
# 
# assume generation time of 10 years (its a perennial plant), 
	# following Cooper et al. Plos ONe 2010 
	# Genetic Variation at Nuclear loci fails to distinguish group is about 3 MY,
		# phy height is 3. So each unit = 1,000,000 years or thus 100,000 generations
# TreeYears=100000
# timeStep<-1/TreeYears
# totalTreeLength=TreeYears*sum(phy$edge.length) #how many generations are represented
# number of expected polinator shifts based on parsimony is 7:
# parsimonyShifts=7
# pollinatorShiftRate=parsimonyShifts/totalTreeLength

# aquilegia regimes - pollinator syndromes
aquilegiaPollinators<-aquilegiaTrait[,14]
# regimes coded 0, 1, 2
	# 0 is bumble-bee, 1 is humming-bird, 2 is hawkmoth
# this probably won't be used directly?
# could use for post-analysis comparisons? Hmm
	

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
# make multiPhylo	
class(ideaTrees)<-"multiPhylo"
# compress tip labels? No, I don't think that works for trees of different sizes
	#	trees<-.compressTipLabel(trees)
		


		
		
		
		
		
		
		
		
		
		
		


# does this analysis-run depend on previous runs or not?
if(dependentPrevRun){

}else{
	
	}




#treeSet

# if the treeSet is "Ideal-Simulated"
# then the number of treeTypes and nTipNumbers is 3, other 1
nTreeTypes<-nTipNumbers<-1
#
if(treeSet == "empirical_anolis_tree"){
	treeList<-
	}
#
if(treeSet == "empirical_Aquilegia_tree"){
	treeList<-
	}
#
if(treeSet=="Ideal_Simulated"){
	treeList<-ideaTrees
	nTreeTypes<-nTipNumbers<-3
	}
#
# empiricalTraitData
#
# nTraitSets is 1 unless empiricalTraitData is "SIMULATED" in which case it is nSimTrait
nTraitSets<-1
if(empiricalTraitData == "Anolis_Size_Data"){
	
	}
#
if(empiricalTraitData == "Aquilegia_Nectar_Spur_Data"){
	
	}
#
if(empiricalTraitData == "SIMULATED"){
	# nTraitSets is 1 unless empiricalTraitData is "SIMULATED" in which case it is nSimTrait
	nTraitSets<-nSimTrait
	# simTrait.Intrinsic
	#
	if(is.na(simTrait.Intrinsic)){
		stop("The intrinsic model for a simulated trait dataset is given as NA")
	}else{
		# ANOLIS BASED MODELS
		if(simTrait.Intrinsic == "An_Emp_BrownMotion"){
			
			}	
		#
		if(simTrait.Intrinsic == "An_Emp_Disp"){
			
			}	
		#
		if(simTrait.Intrinsic == "An_Emp_DispBound"){
			
			}	
		#
		if(simTrait.Intrinsic == "An_Emp_Bound"){
			
			}	
		#
		if(simTrait.Intrinsic == "An_Emp_Bound_BoundByStartingState"){
			
			}	
		#
		if(simTrait.Intrinsic == "An_Emp_Bound_BoundByMinValue"){
			
			}	
		#
		if(simTrait.Intrinsic == "An_Emp_Bound_BoundOneRangeAway"){
			
			}	
		#
		if(simTrait.Intrinsic == "An_Emp_TimeReg"){
			
			}	
		#
		if(simTrait.Intrinsic == "Aq_Emp_3Opt2Bound"){
			
			}	
		#
		if(simTrait.Intrinsic == "Aq_Emp_BrownMotion"){
			
			}	
		#
				
			
		}
	#
	# simTrait.Extrinsic
	#
	if(is.na(simTrait.Extrinsic)){
		stop("The extrinsic model for a simulated trait dataset is given as NA")
	}else{
		
		if(simTrait.Extrinsic == "Null"){
			
			}
		#
		if(simTrait.Extrinsic == "An_Emp_Bound"){
			
			}
		#
		if(simTrait.Extrinsic == "An_Emp_Disp"){
			
			}
		#
		if(simTrait.Extrinsic == "An_Emp_DispBound"){
			
			}
		#
	
		
		}

	
	
	
	}
#
#







doRun.Intrinsic

Pure_BM
BM_w/_LowerBound
3Opt2Bound
Time_AutoRegressive_Model

doRun.Extrinsic

Null
Displacement

prior
standard_(uniform)
rexp_with_mean_NOT_at_true_sigmasq






# calculate the number of doRun statements for this analysis-run
# product of treeTypes and nTipNumbers and nSimTrait
nDoRun <- treeTypes * nTipNumbers * nTraitSets
# should be one 1, 10 or 90... probably
		
		
		
	
	
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





# get simulation run table
simRunTable<-read.csv(header=TRUE,
	stringsAsFactors=FALSE,
	file="simulation_sets_07-20-18.csv")


# remove ID, remove comments
simRunTable<-simRunTable[,-c(1,ncol(simRunTable))]

nRepeatTraitSim<-10			
simNtip<-"c(8, 16, 64)"
idealTreeTypes<-'c("Ideal-Balanced", "Ideal-Pectinate", "Ideal-Star")'

#####################################################

# convert the table
simRunTable$nTipNumbers[simRunTable$nTipNumbers==1]<-NA
simRunTable$nTipNumbers[simRunTable$nTipNumbers==3]<-simNtip
colnames(simRunTable)[colnames(simRunTable)=="nTipNumbers"]<-"nTipSets"
simRunTable$nTreeTypes[simRunTable$nTreeTypes==1]<-NA
simRunTable$nTreeTypes[simRunTable$nTreeTypes==3]<-idealTreeTypes
colnames(simRunTable)[colnames(simRunTable)=="nTreeTypes"]<-"idealTreeSets"

simRunTable$nSimTrait[simRunTable$nSimTrait==1]<-NA

#simRunTable<-simRunTable[,!(colnames(simRunTable)=="nTree")]
simRunTable<-simRunTable[,!(colnames(simRunTable)=="nDoRun")]


	
	

	
	
