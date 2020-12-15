#------------------------------------------------------------------------------------------
# Analyses script of the simulation for the manuscript
#------------------------------------------------------------------------------------------

#Load necessary libraries

library(ape)
library(geiger)
library(phytools)
library(OUwie)
library(dplyr)
library(btw)

setwd("~/ownCloud/Documents/CulturalPhylogenetics/Code/PamaNyungan/")

clade_results<-read.csv("SimulateCulture_MultipleSelectionResults_clade.csv")

ecology_results<-read.csv("SimulateCulture_MultipleSelectionResults_ecology.csv")


drift_variants<-c("drift_two_slow","drift_two_medium","drift_two_fast","drift_four_slow","drift_four_medium","drift_four_fast")
tree_variants <- c("PamaNyungantree","Grafentree","Onetree","Earlytree","Latetree")
selection_variants<-c("selection_two_twofold","selection_two_fourfold","selection_two_eightfold","selection_four_twofold","selection_four_fourfold","selection_four_eightfold", "selection_four_pathway_slow", "selection_four_pathway_fast","selection_four_pathway_fromone")

#LR= 2[log-likelihood(complex model) – log-likelihood(simple model)]
#for drift: compare drift to pahtwaystoone

#dchisq(value,parameters)


clade_results$X.1<-c(1:nrow(clade_results))
ecology_results$X.1<-c(1:nrow(ecology_results))

clade_results[filter(clade_results,(BayestraitsDependent_CladeA-BayestraitsIndependent_CladeA)>50)$X.1,]$BayestraitsDependent_CladeA<-NA
clade_results[filter(clade_results,(BayestraitsIndependent_CladeA-BayestraitsDependent_CladeA)>50)$X.1,]$BayestraitsIndependent_CladeA<-NA
clade_results[filter(clade_results,(BayestraitsIndependent_Forest-BayestraitsDependent_Forest)>50)$X.1,]$BayestraitsIndependent_Forest<-NA
clade_results[filter(clade_results,(BayestraitsDependent_Forest-BayestraitsIndependent_Forest)>50)$X.1,]$BayestraitsDependent_Forest<-NA

ecology_results[filter(ecology_results,(BayestraitsDependent_CladeA-BayestraitsIndependent_CladeA)>50)$X.1,]$BayestraitsDependent_CladeA<-NA
ecology_results[filter(ecology_results,(BayestraitsIndependent_CladeA-BayestraitsDependent_CladeA)>50)$X.1,]$BayestraitsIndependent_CladeA<-NA
ecology_results[filter(ecology_results,(BayestraitsIndependent_Forest-BayestraitsDependent_Forest)>50)$X.1,]$BayestraitsIndependent_Forest<-NA
ecology_results[filter(ecology_results,(BayestraitsDependent_Forest-BayestraitsIndependent_Forest)>50)$X.1,]$BayestraitsDependent_Forest<-NA





#dependent model is the complex one as it has more parameters, Log BF = 2(log marginal likelihood complex model – log marginal likelihood simple model)
clade_results<-mutate(clade_results,BT_dependent_CladeA_better = 2*(clade_results[,]$BayestraitsDependent_CladeA-clade_results[,]$BayestraitsIndependent_CladeA) )
clade_results<-mutate(clade_results,BT_dependent_Forest_better = 2*(clade_results[,]$BayestraitsDependent_Forest-clade_results[,]$BayestraitsIndependent_Forest) )


ecology_results<-mutate(ecology_results,BT_dependent_CladeA_better = 2*(ecology_results[,]$BayestraitsDependent_CladeA-ecology_results[,]$BayestraitsIndependent_CladeA) )
ecology_results<-mutate(ecology_results,BT_dependent_Forest_better = 2*(ecology_results[,]$BayestraitsDependent_Forest-ecology_results[,]$BayestraitsIndependent_Forest) )


#Sample size listing

#Create subsets of trees according to the major language groups
# cladeA  # 102 societies "NorthernAmerind"
# cladeB  # 28 societies "NaDene"
# cladeC  # 42 societies "CentralAmerind"
# cladeD  # 36 societies "Penutian"
# cladeE  # 27 societies "Hokan"
# cladeF  # 39 societies "Almosan"

# Three quarter sample # 129 societies
# Half sample # 86 societies
# State dependent loss of samples - 50% from population of least frequent variant # ~ 167 societies
# State dependent loss of samples - 50% from population of most frequent variant # ~108 societies
# Loss according to ecoregion (50% of those living in forests) # 89 societies live in forest, 45 missing # 127 societies





clade_results<-mutate(clade_results,chisq_BT_CladeAdependent_better=(1-pchisq(clade_results[,]$BT_dependent_CladeA_better, 4 )))
clade_results<-mutate(clade_results,chisq_BT_Forestdependent_better=(1-pchisq(clade_results[,]$BT_dependent_Forest_better, 4 )))

ecology_results<-mutate(ecology_results,chisq_BT_CladeAdependent_better=(1-pchisq(ecology_results[,]$BT_dependent_CladeA_better, 4 )))
ecology_results<-mutate(ecology_results,chisq_BT_Forestdependent_better=(1-pchisq(ecology_results[,]$BT_dependent_Forest_better, 4 )))


clade_results$BT_CladeAdependent_wronglysupported<-as.integer(clade_results$chisq_BT_CladeAdependent_better<0.05)
clade_results$BT_Forestdependent_wronglysupported<-as.integer(clade_results$chisq_BT_Forestdependent_better<0.05)

ecology_results$BT_CladeAdependent_wronglysupported<-as.integer(ecology_results$chisq_BT_CladeAdependent_better<0.05)
ecology_results$BT_Forestdependent_wronglysupported<-as.integer(ecology_results$chisq_BT_Forestdependent_better<0.05)





#write new files to csv
write.csv(clade_results,file="CulturalPhylogenetics_SimulationResults_MultistateClade_ForAnalysis.csv")
write.csv(ecology_results,file="CulturalPhylogenetics_SimulationResults_MultistateEcology_ForAnalysis.csv")









clade_results %>% filter(chisq_BT_CladeAdependent_better>0.05) %>% group_by(Tree) %>% summarise(n_distinct(X))
clade_results %>% filter(chisq_BT_CladeAdependent_better<0.05) %>% group_by(Tree) %>% summarise(n_distinct(X))

clade_results %>% filter(chisq_BT_CladeAdependent_better<0.05) %>% group_by(Sample) %>% summarise(n_distinct(X))
clade_results %>% filter(chisq_BT_CladeAdependent_better>0.05) %>% group_by(Sample) %>% summarise(n_distinct(X))

clade_results %>% filter(chisq_BT_CladeAdependent_better>0.05) %>% group_by(NumberOfVariantsInModel) %>% summarise(n_distinct(X))
clade_results %>% filter(chisq_BT_CladeAdependent_better<0.05) %>% group_by(NumberOfVariantsInModel) %>% summarise(n_distinct(X))

clade_results %>% filter(chisq_BT_CladeAdependent_better>0.05) %>% group_by(RateOfChange) %>% summarise(n_distinct(X))
clade_results %>% filter(chisq_BT_CladeAdependent_better<0.05) %>% group_by(RateOfChange) %>% summarise(n_distinct(X))



ecology_results %>% filter(chisq_BT_Forestdependent_better>0.05) %>% group_by(Tree) %>% summarise(n_distinct(X))
ecology_results %>% filter(chisq_BT_Forestdependent_better<0.05) %>% group_by(Tree) %>% summarise(n_distinct(X))

ecology_results %>% filter(chisq_BT_Forestdependent_better<0.05) %>% group_by(Sample) %>% summarise(n_distinct(X))
ecology_results %>% filter(chisq_BT_Forestdependent_better>0.05) %>% group_by(Sample) %>% summarise(n_distinct(X))

ecology_results %>% filter(chisq_BT_Forestdependent_better>0.05) %>% group_by(NumberOfVariantsInModel) %>% summarise(n_distinct(X))
ecology_results %>% filter(chisq_BT_Forestdependent_better<0.05) %>% group_by(NumberOfVariantsInModel) %>% summarise(n_distinct(X))

ecology_results %>% filter(chisq_BT_Forestdependent_better>0.05) %>% group_by(RateOfChange) %>% summarise(n_distinct(X))
ecology_results %>% filter(chisq_BT_Forestdependent_better<0.05) %>% group_by(RateOfChange) %>% summarise(n_distinct(X))






selection_results %>% filter(horizontaltransmission %in% "No") %>% summarise(n_distinct(X))

selection_results %>% filter(chisq_p_selectioncorrectlysupported>0.05,horizontaltransmission %in% "No") %>% summarise(n_distinct(X))

selection_results %>% filter(chisq_p_selectioncorrectlysupported>0.05,horizontaltransmission %in% "No") %>% group_by(RateOfChange) %>% summarise(n_distinct(X))
selection_results %>% filter(horizontaltransmission %in% "No") %>% group_by(RateOfChange) %>% summarise(n_distinct(X))


selection_results %>% filter(chisq_p_selectioncorrectlysupported>0.05,horizontaltransmission %in% "No") %>% group_by(NumberOfVariantsObserved) %>% summarise(n_distinct(X))
selection_results %>% filter(horizontaltransmission %in% "No") %>% group_by(NumberOfVariantsObserved) %>% summarise(n_distinct(X))


selection_results %>% filter(chisq_p_selectioncorrectlysupported>0.05,horizontaltransmission %in% "Yes") %>% group_by(NumberOfVariantsObserved) %>% summarise(n_distinct(X))
selection_results %>% filter(horizontaltransmission %in% "Yes") %>% group_by(NumberOfVariantsObserved) %>% summarise(n_distinct(X))

selection_results %>% filter(chisq_p_selectioncorrectlysupported>0.05,horizontaltransmission %in% "Yes") %>% group_by(RateOfChange) %>% summarise(n_distinct(X))
selection_results %>% filter(horizontaltransmission %in% "Yes") %>% group_by(RateOfChange) %>% summarise(n_distinct(X))


selection_results %>% filter(chisq_p_selectioncorrectlysupported>0.05,horizontaltransmission %in% "No") %>% group_by(Tree) %>% summarise(n_distinct(X))
selection_results %>% filter(horizontaltransmission %in% "No") %>% group_by(Tree) %>% summarise(n_distinct(X))


selection_results %>% filter(chisq_p_selectioncorrectlysupported>0.05,horizontaltransmission %in% "Yes") %>% group_by(Tree) %>% summarise(n_distinct(X))
selection_results %>% filter(horizontaltransmission %in% "Yes") %>% group_by(Tree) %>% summarise(n_distinct(X))


selection_results %>% filter(chisq_BT_CladeAdependent_better<0.05,horizontaltransmission %in% "No") %>% group_by(Tree) %>% summarise(n_distinct(X))
selection_results %>% filter(chisq_BT_CladeAdependent_better>0.05,horizontaltransmission %in% "No") %>% group_by(Tree) %>% summarise(n_distinct(X))


selection_results %>% filter(chisq_BT_CladeAdependent_better<0.05,horizontaltransmission %in% "Yes") %>% group_by(Tree) %>% summarise(n_distinct(X))
selection_results %>% filter(chisq_BT_CladeAdependent_better>0.05,horizontaltransmission %in% "Yes") %>% group_by(Tree) %>% summarise(n_distinct(X))


selection_results %>% filter(chisq_BT_Forestdependent_better<0.05,horizontaltransmission %in% "No") %>% group_by(Tree) %>% summarise(n_distinct(X))
selection_results %>% filter(chisq_BT_Forestdependent_better>0.05,horizontaltransmission %in% "No") %>% group_by(Tree) %>% summarise(n_distinct(X))


drift_results %>% filter(chisq_BT_CladeAdependent_better<0.05,horizontaltransmission %in% "Yes") %>% group_by(Tree) %>% summarise(n_distinct(X))
drift_results %>% filter(chisq_BT_CladeAdependent_better>0.05,horizontaltransmission %in% "Yes") %>% group_by(Tree) %>% summarise(n_distinct(X))


drift_results %>% filter(chisq_BT_Forestdependent_better<0.05,horizontaltransmission %in% "No") %>% group_by(Tree) %>% summarise(n_distinct(X))
drift_results %>% filter(chisq_BT_Forestdependent_better>0.05,horizontaltransmission %in% "No") %>% group_by(Tree) %>% summarise(n_distinct(X))



selection_results %>% filter(chisq_BT_CladeAdependent_better<0.05) %>% group_by(Sample) %>% summarise(n_distinct(X))
selection_results %>% filter(chisq_BT_CladeAdependent_better>0.05) %>% group_by(Sample) %>% summarise(n_distinct(X))

drift_results %>% filter(chisq_BT_CladeAdependent_better<0.05) %>% group_by(Sample) %>% summarise(n_distinct(X))
drift_results %>% filter(chisq_BT_CladeAdependent_better>0.05) %>% group_by(Sample) %>% summarise(n_distinct(X))





plot(samplesize_analysis[1:12,]$proportion_drift_selectionwronglysupported~samplesize_analysis[1:12,]$samplesize,xlim=c(0,200),ylim=c(0,1))
points(samplesize_analysis$proportion_selection_selectionwronglynotsupported~samplesize_analysis$samplesize,col="red")
points(samplesize_analysis$drift_diversityrepresented~samplesize_analysis$samplesize,col="black",pch=8)
points(samplesize_analysis$selection_diversityrepresented~samplesize_analysis$samplesize,col="red",pch=8)




#outcome variables: 
#drift: proportion wrongly assigned as selection (~20%) - false positives
#drift: proportion of total variants recovered
#drift: Clade A Bayestraits wrongly supported
#drift: Forest Bayestraits wrongly supported
#drift: Bayestraits Multistate dependent wrongly supported

#selection: proportion wrongly not detected as selection (~50%) - false negatives
#selection: proportion of total variants recovered
#selection: Clade A Bayestraits wrongly supported
#selection: Forest Bayestraits wrongly supported
#selection: Bayestraits Multistate dependent wrongly supported

#predictor variables:
#sample size
#rate of change
#number of variants in model
#horizontal transmission
#tree shape
#for selection: to one or pathway



drift_results_filtered<-drift_results[ complete.cases(drift_results$chisq_p_selectionwronglysupported) , ]


dat_list <- list(
  chisq_p_selectionwronglysupported = drift_results_filtered$chisq_p_selectionwronglysupported,
  samplesize = standardize(drift_results_filtered$samplesize),
  horizontaltransmission = as.integer(drift_results_filtered$horizontaltransmission)
)

m8.1 <- quap( alist(
  chisq_p_selectionwronglysupported ~ dnorm( mu , sigma ) , 
  mu <- a + b*samplesize , 
  a ~ dnorm( 1 , 1 ) ,
  b ~ dnorm( 0 , 1 ) ,
  sigma ~ dexp( 1 ) ) , data=dat_list )
precis(m8.1)


set.seed(7)
prior <- extract.prior( m8.1 )
# set up the plot dimensions
plot( NULL , xlim=c(0,1) , ylim=c(0.5,1.5) ,
      xlab="ruggedness" , ylab="log GDP" ) abline( h=min(dd$log_gdp_std) , lty=2 ) abline( h=max(dd$log_gdp_std) , lty=2 )
# draw 50 lines from the prior
rugged_seq <- seq( from=-0.1 , to=1.1 , length.out=30 )
mu <- link( m8.1 , post=prior , data=data.frame(rugged_std=rugged_seq) ) for ( i in 1:50 ) lines( rugged_seq , mu[i,] , col=col.alpha("black",0.3) )


Latitudedistance<-as.matrix(dist(WNAIlocations$Latitude))
Longitudedistance<-as.matrix(dist(WNAIlocations$Longitude))
Totaldistance<-Latitudedistance+Longitudedistance
diag(Totaldistance)<-0
colnames(Totaldistance)<-row.names(WNAIlocations)
row.names(Totaldistance)<-row.names(WNAIlocations)
distancetree = nj(Totaldistance)
plot.phylo(distancetree,edge.width=3,label.offset=1, cex=0.4)



plot.phylo(Grafentree,type="fan",edge.color = (simmapattemptEcology$mapped.edge[,1]>0.5)+1,show.tip.label = F,edge.width=3)

plot.phylo(Onetree,edge.width=3,label.offset=1, cex=0.4,direction="leftwards")

