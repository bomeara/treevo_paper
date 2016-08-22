library(OUwie)
source("/Users/bomeara/Documents/MyDocuments/Active/OUwie/pkg/R/OUwie.fixed.R")
source("/Users/bomeara/Documents/MyDocuments/Active/OUwie/pkg/R/varcov.ou.R")
source("/Users/bomeara/Documents/MyDocuments/Active/OUwie/pkg/R/weight.mat.R")
library(TreEvo)
TreEvo.files<-system("ls -1 /Users/bomeara/Documents/MyDocuments/Active/TreEvo/pkg/R/*.R",intern=TRUE)
for (i in sequence(length(TreEvo.files))) {
  source(TreEvo.files[i]) 
}
abcTolerance<-1

likelihoodBM.fixed<-function(phy,x,rate,root.state) {
  k<-2
  
  vcv<-vcv.phylo(phy)
  
  
    vv<-exp(x)*vcv
    mu<-root.state
    mu<-rep(mu, n)
    return(-dmvnorm(y, mu, vv, log=T))

  
}

setwd("/Users/bomeara/Documents/MyDocuments/Ongoing/Talks/2012vii07_Evolution2012_TreEvo/Simulation_BMSimple2_runL")
load("ntax30_bm_rej.intialsettings.Rsave")
load("ntax30_bm_rej.trueFreeValuesANDSummaryValues.Rsave")
trueFreeValuesMatrix<-trueFreeValuesANDSummaryValues[,1:2] #NOTE I AM TAKING A SUBSET AS I AM BORED
summaryValuesMatrix<-trueFreeValuesANDSummaryValues[,c(-1,-2)]
colnames(traits)<-"X"
OUwie.phy<-phy
OUwie.phy$node.label<-1+rbinom(Nnode(phy),1,0.5)
OUwie.traits<-data.frame(Genus_species=rownames(traits),Regime=rep(1,30),X=traits)
a<-OUwie(OUwie.phy, OUwie.traits,model="BM1")
#BM1 won't use root states, so have to use OU equivalent with VERY weak alpha
b<-OUwie.fixed(OUwie.phy, OUwie.traits,model="OU1",sigma.sq=15.4785, root.station=FALSE, alpha=0.00000001, theta=rep(8.749418,2))

abcTolerance<-1
res<-PLSRejection(summaryValuesMatrix, trueFreeValuesMatrix, phy, traits, abcTolerance)

#Note that we have uniform priors for both params

#ouwie scales its tree to height 1
numberofsteps<-floor(splits[1, 1]/timeStep)
ScaleRate<-function(discreteRate, numberofsteps) {
  sigma.sq<-(discreteRate^2)*numberofsteps
}

results<-cbind(res[[1]],rep(NA,dim(res[[1]])[1]))

save(results,file="results_with_distances.RSave",compress=TRUE)


for (i in sequence(dim(results)[1])) {
  results[i,dim(results)[2]]<-OUwie.fixed(OUwie.phy, OUwie.traits,model="OU1",sigma.sq=ScaleRate(results$param2[i],numberofsteps), root.station=FALSE, alpha=0.00000001, theta=rep(results$param1[i],2))$loglik
  if (i%%100==0) {
  	print(i)
  	save(results,file="results_with_likelihood_and_distances.RSave",compress=TRUE)
  }
  #print(results[i,])
}


colnames(results)[9]<-"loglik"
#par(mfrow=c(2,2))
#results.pruned<-results[which(results$loglik>quantile(results$loglik,probs=0.5)),] #get rid of weird outliers
#plot(results.pruned$param1, results.pruned$loglik)


#plot(results.pruned$param2, 1-(results.pruned$loglik)/min((results.pruned$loglik)),col="red",ylim=c(0,1))

#my.lines<-density((results.pruned[which(results.pruned$distance<quantile(results.pruned$distance,1)),]$param2))
#lines(my.lines$x,my.lines$y/max(my.lines$y))

#my.lines<-density((results.pruned[which(results.pruned$distance<quantile(results.pruned$distance,.1)),]$param2))
#lines(my.lines$x,my.lines$y/max(my.lines$y))

#my.lines<-density((results.pruned[which(results.pruned$distance<quantile(results.pruned$distance,.01)),]$param2))
#lines(my.lines$x,my.lines$y/max(my.lines$y))



#plot(density((results.pruned[which(results.pruned$distance<quantile(results.pruned$distance,0.1)),]$param1)),col="red")
#(hist((results.pruned[which(results.pruned$distance<quantile(results.pruned$distance,1)),]$param2)))
#(hist((results.pruned[which(results.pruned$distance<quantile(results.pruned$distance,0.1)),]$param2)))
#(hist((results.pruned[which(results.pruned$distance<quantile(results.pruned$distance,0.01)),]$param2)))
#(hist((results.pruned[which(results.pruned$distance<quantile(results.pruned$distance,0.001)),]$param2)))



#plot(results.pruned$distance, -results.pruned$loglik,log="x")
save(results,file="results_with_likelihood_and_distances.RSave",compress=TRUE)
best.param2<-results.pruned[which.max(results.pruned$loglik),]$param2

par(mfrow=c(2,2))
plot(results$param1, results$loglik-max(results$loglik),ylim=c(-10,0),col="grey",pch=".")
lines(x=range(results[which(results$loglik-max(results$loglik)>-2),]$param1),y=rep(-2,2),col="red")
cutoff=0.001
lines(x=range(results[which(results$distance<quantile(results$distance,cutoff)) ,]$param1),y=rep(-cutoff,2))

plot(results$param2, results$loglik-max(results$loglik),ylim=c(-10,0),col="grey",pch=".")
lines(x=range(results[which(results$loglik-max(results$loglik)>-2),]$param2),y=rep(-2,2),col="red")
cutoff=0.001
lines(x=range(results[which(results$distance<quantile(results$distance,cutoff)) ,]$param2),y=rep(-cutoff,2))


cutoffs<-10^(0:-6)
plot(range(cutoffs),range(results$param1),type="n",bty="n",log="x") 
for(i in sequence(length(cutoffs))) {
  cutoff<-cutoffs[i]
  lines(y=range(results[which(results$distance<quantile(results$distance,cutoff)) ,]$param1),x=rep(cutoffs[i],2))
  
}
lines(y=rep(max(results[which(results$loglik-max(results$loglik)>-2),]$param1),2),x=range(cutoffs),col="red")
lines(y=rep(min(results[which(results$loglik-max(results$loglik)>-2),]$param1),2),x=range(cutoffs),col="red")

plot(range(cutoffs),range(results$param2),type="n",bty="n",log="x") 
for(i in sequence(length(cutoffs))) {
  cutoff<-cutoffs[i]
  lines(y=range(results[which(results$distance<quantile(results$distance,cutoff)) ,]$param2),x=rep(cutoffs[i],2))
  
}
lines(y=rep(max(results[which(results$loglik-max(results$loglik)>-2),]$param2),2),x=range(cutoffs),col="red")
lines(y=rep(min(results[which(results$loglik-max(results$loglik)>-2),]$param2),2),x=range(cutoffs),col="red")

