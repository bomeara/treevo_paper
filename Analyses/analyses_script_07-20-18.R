###################################################
# Individual Empirical Analyses and Simulations
##################################################

setwd("d://dave//workspace//treevo_paper//")

library(ape)
library(TreEvo)

######################################
# get empirical data

# obtain anolis tree - from Poe et al. 2017 (SystBiol)
	# their time-tree
anolisTree<-read.tree(file="datasets//anolis_PoeEtAl2018_datedMCCw0.5burnin.tre")

# obtain anolis trait data - 
	# Snout-Vent body-size data from Poe et al. 2018 (AmNat)
anolisTrait<-read.table("datasets//anolis_lntraits_matched_tabdelim_07-24-18.txt",
	header=TRUE,row.names=1)
anolisSize<-anolisTrait[,1]

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
aquilegiaSpurLength<-ln(aquilegiaSpurLength)

# aquilegia regimes - pollinator syndromes
aquilegiaPollinators<-aquilegiaTrait[,14]
# regimes coded 0, 1, 2
	# 0 is bumble-bee, 1 is humming-bird, 2 is hawkmoth
	
	
######	
oldCode<-'	
setwd("/Users/bomeara/Documents/MyDocuments/Ongoing/Talks/2012vii07_Evolution2012_TreEvo/Aquilegia")
library(TreEvo)
phy<-read.tree("Aquilegia.phy")[[1]]
data<-read.table("Aquilegia.txt",header=F,sep="\t")
char<-log(data[,3]) #log spur length
abcMethod<-"rejection"
abcTolerance<-0.01
library(doMC, quietly=T)
library(foreach, quietly=T)
library(parallel)
registerDoMC()

pollinatorShiftIntrinsic <-function(params, states, timefrompresent) {
  #params[1] is sd for normal process, 
  #params[2] is sd for jump size of shift, 
  #params[3] is mean of shift (mean of normal process is zero), 
  #params[4] is per generation probability of shift
  background.process<-TRUE
  if (runif(1,0,1)<params[4]) {
    background.process<-FALSE 
  }
  if (background.process) {
    newdisplacement<-rnorm(n=length(states),mean=0,sd=params[1])
    return(newdisplacement)
  }
  else {
    newdisplacement<-rnorm(n=length(states),mean=params[3],sd=params[2])
    return(newdisplacement)
  }
}

#assume generation time of 10 years (its a perennial plant), following Cooper et al. Plos ONe 2010 Genetic Variation at Nuclear loci fails to distinguish
#group is about 3 MY, phy height is 3. So each unit = 1,000,000 years or thus 100,000 generations
TreeYears=100000
timeStep<-1/TreeYears

totalTreeLength=TreeYears*sum(phy$edge.length) #how many generations are represented
parsimonyShifts=7
pollinatorShiftRate=parsimonyShifts/totalTreeLength

maxPatristicDistance<-max(branching.times(phy))*TreeYears
maxCharDistance<-max(char)-min(char)
predictDistance<-function(x,maxPatristicDistance=maxPatristicDistance,nreps=20) {
  return(mean(replicate(nreps,abs(sum(rnorm(round(maxPatristicDistance),mean=0,sd=x))))))
}
distanceDistance<-function(x,maxPatristicDistance=maxPatristicDistance,maxCharDistance,nreps=20) {
   return(abs(maxCharDistance - predictDistance(x,maxPatristicDistance,nreps)))
}


#interval=c(0.0001,.05)
#nreps=200
#results.optimize<-optimize(f=distanceDistance,interval=interval,maxPatristicDistance=maxPatristicDistance,maxCharDistance=maxCharDistance,nreps=nreps)
#sd.bestguess<-results.optimize$minimum
sd.bestguess<-0.008 #it is estimated as this from above procedure. This way, all sims use same priors rather than reoptimizing

#x<-seq(from=min(interval),to=max(interval),length.out=10)
#plot(x,sapply(x,distanceDistance,maxPatristicDistance=maxPatristicDistance,maxCharDistance=maxCharDistance,nreps=nreps))

intrinsicFn=pollinatorShiftIntrinsic
extrinsicFn=nullExtrinsic
startingPriorsFns="fixed"
startingPriorsValues=matrix(rep(min(char),2),nrow=2,byrow=FALSE) #assume that the min value is the root state
intrinsicPriorsFns=c("exponential","fixed","uniform","uniform") #do fixed for param2 b/c dont have a lot of data to use
intrinsicPriorsValues=matrix(c(
  rep(1/sd.bestguess , 2),
  0, 0 ,
  -maxCharDistance, maxCharDistance,
  0, 0.00001
  ), nrow=2, byrow=FALSE)
extrinsicPriorsFns="fixed"
extrinsicPriorsValues=matrix(c(0, 0), nrow=2, byrow=FALSE)
	
'	

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
	
#############################################################################
# analyses


#############################################
   # runLabel = An_Emp-DispBound 
   # treeSet = empirical-anolis_tree 
   # empiricalTraitData = Anolis_Size_Data 
   # doRun.Intrinsic = BM_w/_LowerBound 
   # doRun.Extrinsic = Displacement 
   # prior = standard_(uniform) 


intrinsicModel<-
extrinsicModel<-
   
   



#############################################
   # runLabel = An_Emp-Bound 
   # treeSet = empirical-anolis_tree 
   # empiricalTraitData = Anolis_Size_Data 
   # doRun.Intrinsic = BM_w/_LowerBound 
   # doRun.Extrinsic = Null 
   # prior = standard_(uniform) 






#############################################
   # runLabel = An_SimDispBound-DispBound 
   # treeSet = empirical-anolis_tree 
   # simTrait.Intrinsic = An_Emp-DispBound 
   # simTrait.Extrinsic = An_Emp-DispBound 
   # doRun.Intrinsic = BM_w/_LowerBound 
   # doRun.Extrinsic = Displacement 
   # prior = standard_(uniform) 
   # nSimTrait = 10 






#############################################
   # runLabel = An_SimBound-DispBound 
   # treeSet = empirical-anolis_tree 
   # simTrait.Intrinsic = An_Emp-Bound 
   # simTrait.Extrinsic = Null 
   # doRun.Intrinsic = Pure_BM 
   # doRun.Extrinsic = Displacement 
   # prior = standard_(uniform) 
   # nSimTrait = 10 






#############################################
   # runLabel = An_Emp-BrownMotion 
   # treeSet = empirical-anolis_tree 
   # empiricalTraitData = Anolis_Size_Data 
   # doRun.Intrinsic = Pure_BM 
   # doRun.Extrinsic = Null 
   # prior = standard_(uniform) 






#############################################
   # runLabel = An_Emp-Disp 
   # treeSet = empirical-anolis_tree 
   # empiricalTraitData = Anolis_Size_Data 
   # doRun.Intrinsic = Pure_BM 
   # doRun.Extrinsic = Displacement 
   # prior = standard_(uniform) 






#############################################
   # runLabel = An_Emp-TimeReg 
   # treeSet = empirical-anolis_tree 
   # empiricalTraitData = Anolis_Size_Data 
   # doRun.Intrinsic = Pure_BM 
   # doRun.Extrinsic = Null 
   # prior = standard_(uniform) 






#############################################
   # runLabel = Aq_Emp-3Opt 
   # treeSet = empirical-Aquilegia_tree 
   # empiricalTraitData = Aquilegia_Nectar_Spur_Data 
   # doRun.Intrinsic = 3-Optima 
   # doRun.Extrinsic = Null 
   # prior = standard_(uniform) 


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
   




#############################################
   # runLabel = Aq_Emp-BrownMotion 
   # treeSet = empirical-Aquilegia_tree 
   # empiricalTraitData = Aquilegia_Nectar_Spur_Data 
   # doRun.Intrinsic = Pure_BM 
   # doRun.Extrinsic = Null 
   # prior = standard_(uniform) 






#############################################
   # runLabel = Aq_Sim3Opt-3Opt 
   # treeSet = empirical-Aquilegia_tree 
   # simTrait.Intrinsic = Aq_Emp-3Opt 
   # simTrait.Extrinsic = Null 
   # doRun.Intrinsic = 3-Optima 
   # doRun.Extrinsic = Null 
   # prior = standard_(uniform) 
   # nSimTrait = 10 






#############################################
   # runLabel = Aq_SimBM-3Opt 
   # treeSet = empirical-Aquilegia_tree 
   # simTrait.Intrinsic = Aq_Emp-BrownMotion 
   # simTrait.Extrinsic = Null 
   # doRun.Intrinsic = 3-Optima 
   # doRun.Extrinsic = Null 
   # prior = standard_(uniform) 
   # nSimTrait = 10 






#############################################
   # runLabel = Ideal_SimBM-BM 
   # treeSet = Ideal-Simulated 
   # simTrait.Intrinsic = An_Emp-BrownMotion 
   # simTrait.Extrinsic = Null 
   # doRun.Intrinsic = Pure_BM 
   # doRun.Extrinsic = Null 
   # prior = standard_(uniform) 
   # idealTreeSets = c("Ideal-Balanced", "Ideal-Pectinate", "Ideal-Star") 
   # nTipSets = c(8, 16, 64) 
   # nSimTrait = 10 






#############################################
   # runLabel = Ideal_SimBMpriorBiased 
   # treeSet = Ideal-Simulated 
   # simTrait.Intrinsic = An_Emp-BrownMotion 
   # simTrait.Extrinsic = Null 
   # doRun.Intrinsic = Pure_BM 
   # doRun.Extrinsic = Null 
   # prior = rexp_with_mean_*not*_at_true_sigmasq 
   # idealTreeSets = c("Ideal-Balanced", "Ideal-Pectinate", "Ideal-Star") 
   # nTipSets = c(8, 16, 64) 
   # nSimTrait = 10 






#############################################
   # runLabel = Ideal_SimBM-Disp 
   # treeSet = Ideal-Simulated 
   # simTrait.Intrinsic = An_Emp-BrownMotion 
   # simTrait.Extrinsic = Null 
   # doRun.Intrinsic = Pure_BM 
   # doRun.Extrinsic = Displacement 
   # prior = standard_(uniform) 
   # idealTreeSets = c("Ideal-Balanced", "Ideal-Pectinate", "Ideal-Star") 
   # nTipSets = c(8, 16, 64) 
   # nSimTrait = 10 






#############################################
   # runLabel = Ideal_SimDisp-Disp 
   # treeSet = Ideal-Simulated 
   # simTrait.Intrinsic = An_Emp-Disp 
   # simTrait.Extrinsic = An_Emp-Disp 
   # doRun.Intrinsic = Pure_BM 
   # doRun.Extrinsic = Displacement 
   # prior = standard_(uniform) 
   # idealTreeSets = c("Ideal-Balanced", "Ideal-Pectinate", "Ideal-Star") 
   # nTipSets = c(8, 16, 64) 
   # nSimTrait = 10 






#############################################
   # runLabel = Ideal_SimBound-DispBound 
   # treeSet = Ideal-Simulated 
   # simTrait.Intrinsic = An_Emp-Bound 
   # simTrait.Extrinsic = Null 
   # doRun.Intrinsic = BM_w/_LowerBound 
   # doRun.Extrinsic = Displacement 
   # prior = standard_(uniform) 
   # idealTreeSets = c("Ideal-Balanced", "Ideal-Pectinate", "Ideal-Star") 
   # nTipSets = c(8, 16, 64) 
   # nSimTrait = 10 






#############################################
   # runLabel = Ideal_SimDispBoundNear-DispBound 
   # treeSet = Ideal-Simulated 
   # simTrait.Intrinsic = An_Emp-Bound_BoundByStartingState 
   # simTrait.Extrinsic = An_Emp-Bound 
   # doRun.Intrinsic = BM_w/_LowerBound 
   # doRun.Extrinsic = Displacement 
   # prior = standard_(uniform) 
   # idealTreeSets = c("Ideal-Balanced", "Ideal-Pectinate", "Ideal-Star") 
   # nTipSets = c(8, 16, 64) 
   # nSimTrait = 10 






#############################################
   # runLabel = Ideal_SimDispBoundMod-DispBound 
   # treeSet = Ideal-Simulated 
   # simTrait.Intrinsic = An_Emp-Bound_BoundByMinValue 
   # simTrait.Extrinsic = An_Emp-Bound 
   # doRun.Intrinsic = BM_w/_LowerBound 
   # doRun.Extrinsic = Displacement 
   # prior = standard_(uniform) 
   # idealTreeSets = c("Ideal-Balanced", "Ideal-Pectinate", "Ideal-Star") 
   # nTipSets = c(8, 16, 64) 
   # nSimTrait = 10 






#############################################
   # runLabel = Ideal_SimDispBoundFar-DispBound 
   # treeSet = Ideal-Simulated 
   # simTrait.Intrinsic = An_Emp-Bound_BoundOneRangeAway 
   # simTrait.Extrinsic = An_Emp-Bound 
   # doRun.Intrinsic = BM_w/_LowerBound 
   # doRun.Extrinsic = Displacement 
   # prior = standard_(uniform) 
   # idealTreeSets = c("Ideal-Balanced", "Ideal-Pectinate", "Ideal-Star") 
   # nTipSets = c(8, 16, 64) 
   # nSimTrait = 10 






#############################################
   # runLabel = Ideal_SimDispBound-DispBound 
   # treeSet = Ideal-Simulated 
   # simTrait.Intrinsic = An_Emp-DispBound 
   # simTrait.Extrinsic = An_Emp-DispBound 
   # doRun.Intrinsic = BM_w/_LowerBound 
   # doRun.Extrinsic = Displacement 
   # prior = standard_(uniform) 
   # idealTreeSets = c("Ideal-Balanced", "Ideal-Pectinate", "Ideal-Star") 
   # nTipSets = c(8, 16, 64) 
   # nSimTrait = 10 






#############################################
   # runLabel = Ideal_SimBM-TimeReg 
   # treeSet = Ideal-Simulated 
   # simTrait.Intrinsic = An_Emp-BrownMotion 
   # simTrait.Extrinsic = Null 
   # doRun.Intrinsic = Time-AutoRegressive_Model 
   # doRun.Extrinsic = Null 
   # prior = standard_(uniform) 
   # idealTreeSets = c("Ideal-Balanced", "Ideal-Pectinate", "Ideal-Star") 
   # nTipSets = c(8, 16, 64) 
   # nSimTrait = 10 






#############################################
   # runLabel = Ideal_SimTimeReg-TimeReg 
   # treeSet = Ideal-Simulated 
   # simTrait.Intrinsic = An_Emp-TimeReg 
   # simTrait.Extrinsic = Null 
   # doRun.Intrinsic = Time-AutoRegressive_Model 
   # doRun.Extrinsic = Null 
   # prior = standard_(uniform) 
   # idealTreeSets = c("Ideal-Balanced", "Ideal-Pectinate", "Ideal-Star") 
   # nTipSets = c(8, 16, 64) 
   # nSimTrait = 10 






#############################################
   # runLabel = Ideal_SimBM-3Opt 
   # treeSet = Ideal-Simulated 
   # simTrait.Intrinsic = Aq_Emp-BrownMotion 
   # simTrait.Extrinsic = Null 
   # doRun.Intrinsic = 3-Optima 
   # doRun.Extrinsic = Null 
   # prior = standard_(uniform) 
   # idealTreeSets = c("Ideal-Balanced", "Ideal-Pectinate", "Ideal-Star") 
   # nTipSets = c(8, 16, 64) 
   # nSimTrait = 10 






#############################################
   # runLabel = Ideal_Sim3Opt-3Opt 
   # treeSet = Ideal-Simulated 
   # simTrait.Intrinsic = Aq_Emp-3Opt 
   # simTrait.Extrinsic = Null 
   # doRun.Intrinsic = 3-Optima 
   # doRun.Extrinsic = Null 
   # prior = standard_(uniform) 
   # idealTreeSets = c("Ideal-Balanced", "Ideal-Pectinate", "Ideal-Star") 
   # nTipSets = c(8, 16, 64) 
   # nSimTrait = 10 
