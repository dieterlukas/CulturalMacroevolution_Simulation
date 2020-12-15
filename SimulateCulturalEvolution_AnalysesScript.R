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

drift_results<-read.csv("SimulateCulture_DriftResults.csv")


drift_variants<-c("drift_two_slow","drift_two_medium","drift_two_fast","drift_four_slow","drift_four_medium","drift_four_fast")
tree_variants <- c("PamaNyungantree","Grafentree","Onetree","Earlytree","Latetree")
selection_variants<-c("selection_two_twofold","selection_two_fourfold","selection_two_eightfold","selection_four_twofold","selection_four_fourfold","selection_four_eightfold", "selection_four_pathway_slow", "selection_four_pathway_fast","selection_four_pathway_fromone")


#need to make "selection_variant" variable based on information in "NumberOfVariatnsInModel" and "RateOfChange"
#drift_variants<-c("drift_two_slow","drift_two_medium","drift_two_fast","drift_four_slow","drift_four_medium","drift_four_fast")


#need to replace tree numbers with tree_variants <- c("Grafentree","Onetree","Earlytree","Latetree")
drift_results<-drift_results
drift_results<-mutate(drift_results,Tree=tree_variants[Tree])
drift_results$selection_variant<-"NA"
drift_results[filter(drift_results,NumberOfVariantsInModel %in% 2,RateOfChange %in% 1)$X,]$selection_variant<-drift_variants[1]
drift_results[filter(drift_results,NumberOfVariantsInModel %in% 2,RateOfChange %in% 2)$X,]$selection_variant<-drift_variants[2]
drift_results[filter(drift_results,NumberOfVariantsInModel %in% 2,RateOfChange %in% 3)$X,]$selection_variant<-drift_variants[3]
drift_results[filter(drift_results,NumberOfVariantsInModel %in% 4,RateOfChange %in% 1)$X,]$selection_variant<-drift_variants[4]
drift_results[filter(drift_results,NumberOfVariantsInModel %in% 4,RateOfChange %in% 2)$X,]$selection_variant<-drift_variants[5]
drift_results[filter(drift_results,NumberOfVariantsInModel %in% 4,RateOfChange %in% 3)$X,]$selection_variant<-drift_variants[6]




selection_results<-read.csv("SimulateCulture_SelectionResults.csv")


colnames(drift_results)[23]<-"simulation_variant"
colnames(selection_results)[23]<-"simulation_variant"

drift_results$ml_model_supported<-"NA"
drift_results[drift_results$NumberOfVariantsObserved %in% c(2,3,4),]$ml_model_supported<-colnames(drift_results[drift_results$NumberOfVariantsObserved  %in% c(2,3,4),8:12])[apply(drift_results[drift_results$NumberOfVariantsObserved  %in% c(2,3,4),8:12],1,which.max)]

selection_results$ml_model_supported<-"NA"
selection_results[selection_results$NumberOfVariantsObserved %in% c(2,3,4),]$ml_model_supported<-colnames(selection_results[selection_results$NumberOfVariantsObserved  %in% c(2,3,4),8:12])[apply(selection_results[selection_results$NumberOfVariantsObserved  %in% c(2,3,4),8:12],1,which.max)]

#LR= 2[log-likelihood(complex model) – log-likelihood(simple model)]
#for drift: compare drift to pahtwaystoone

#dchisq(value,parameters)

drift_results<-drift_results[,1:24]

#Likelihood ratio test to assess whether unconstrained model is significantly better than equal rates drift model:
drift_results<-mutate(drift_results,chisq_p_selectionwronglysupported_df = 1-pchisq(2*(drift_results[,]$loglik_unconstrained-drift_results[,]$loglik_drift),   drift_results[,]$NumberOfVariantsObserved*(drift_results[,]$NumberOfVariantsObserved-1)/2   ) )


drift_results<-mutate(drift_results,chisq_p_selectionwronglysupported = 1-pchisq(2*(drift_results[,]$loglik_unconstrained-drift_results[,]$loglik_drift),   drift_results[,]$NumberOfVariantsObserved*(drift_results[,]$NumberOfVariantsObserved-1)/2   ) )
drift_results<-mutate(drift_results,chisq_p_selectiontoonestatewronglysupported = 1-pchisq(2*(drift_results[,]$loglik_selectiononestate-drift_results[,]$loglik_drift),   drift_results[,]$NumberOfVariantsObserved-1   ) )
drift_results<-drift_results[,4:36]
drift_results$selectionwronglysupportedsignificant<-as.integer(drift_results$chisq_p_selectionwronglysupported<0.05)
drift_results$selectiontoonewronglysupportedsignificant<-as.integer(drift_results$chisq_p_selectiontoonestatewronglysupported<0.05)

selection_results<-mutate(selection_results,chisq_p_selectionwcorrectlysupported = 1-pchisq(2*(selection_results[,]$loglik_unconstrained-selection_results[,]$loglik_drift),   selection_results[,]$NumberOfVariantsObserved*(selection_results[,]$NumberOfVariantsObserved-1)/2   ) )
selection_results<-mutate(selection_results,chisq_p_selectiontoonestatecorrectlysupported = 1-pchisq(2*(selection_results[,]$loglik_selectiononestate-selection_results[,]$loglik_drift),   selection_results[,]$NumberOfVariantsObserved-1   ) )
selection_results<-selection_results[,4:34]
selection_results$selectionwronglysupportedsignificant<-as.integer(selection_results$chisq_p_selectionwcorrectlysupported<0.05)
selection_results$selectiontoonewronglysupportedsignificant<-as.integer(selection_results$chisq_p_selectiontoonestatecorrectlysupported<0.05)



#Likelihood ratio test to assess whether directed selection model is significantly better than equal rates drift model:
drift_results<-mutate(drift_results,chisq_p_selectiontoonestatewronglysupported_df = 1-pchisq(2*(drift_results[,]$loglik_selectiononestate-drift_results[,]$loglik_drift), drift_results[,]$NumberOfVariantsObserved-1) )

nrow(filter(drift_results,chisq_p_selectiontoonestatewronglysupported_df < 0.05))
nrow(drift_results)

nrow(filter(drift_results,chisq_p_selectiontoonestatewronglysupported_df < 0.05))/nrow(drift_results)

#dependent model is the complex one as it has more parameters, Log BF = 2(log marginal likelihood complex model – log marginal likelihood simple model)
drift_results<-mutate(drift_results,BT_dependent_CladeA_better = 2*(drift_results[,]$BayestraitsDependent_CladeA-drift_results[,]$BayestraitsIndependent_CladeA) )
drift_results<-mutate(drift_results,BT_dependent_Forest_better = 2*(drift_results[,]$BayestraitsDependent_Ecology-drift_results[,]$BayestraitsIndependent_Ecology) )
drift_results<-mutate(drift_results,BT_dependent_Forest_better = 2*(drift_results[,]$BayestraitsMultistate_Free-drift_results[,]$BayestraitsMultistate_Equal) )


nrow(filter(drift_results,BT_dependent_CladeA_better > 2))
nrow(filter(drift_results,BT_dependent_CladeA_better > 5))

nrow(filter(drift_results,BT_dependent_CladeA_better > 2))/nrow(drift_results)

nrow(filter(drift_results,BT_dependent_Forest_better > 2))
nrow(filter(drift_results,BT_dependent_Forest_better > 5))

nrow(filter(drift_results,BT_dependent_Forest_better > 2))/nrow(drift_results)


#--#--#--#--#--#--#--#--#-#--#

#Likelihood ratio test to assess whether unconstrained model is significantly better than equal rates drift model:
selection_results<-mutate(selection_results,chisq_p_selectionwcorrectlysupported = 1-pchisq(2*(selection_results[,]$loglik_unconstrained-selection_results[,]$loglik_drift),1) )

#Likelihood ratio test to assess whether directed selection model is significantly better than equal rates drift model:
selection_results<-mutate(selection_results,chisq_p_selectiontoonestatecorrectlysupported = 1-pchisq(2*(selection_results[,]$loglik_selectiononestate-selection_results[,]$loglik_drift),1) )

nrow(filter(selection_results,chisq_p_selectionwcorrectlysupported < 0.05))
nrow(selection_results)

nrow(filter(selection_results,chisq_p_selectionwcorrectlysupported < 0.05))/nrow(selection_results)

#dependent model is the complex one as it has more parameters, Log BF = 2(log marginal likelihood complex model – log marginal likelihood simple model)
selection_results<-mutate(selection_results,BT_dependent_CladeA_better = 2*(selection_results[,]$BayestraitsDependent_CladeA-selection_results[,]$BayestraitsIndependent_CladeA) )
selection_results<-mutate(selection_results,BT_dependent_Forest_better = 2*(selection_results[,]$BayestraitsDependent_Forest-selection_results[,]$BayestraitsIndependent_Forest) )

nrow(filter(selection_results,BT_dependent_CladeA_better > 2))
nrow(filter(selection_results,BT_dependent_CladeA_better > 5))

nrow(filter(selection_results,BT_dependent_CladeA_better > 2))/nrow(selection_results)

nrow(filter(selection_results,BT_dependent_Forest_better > 2))
nrow(filter(selection_results,BT_dependent_Forest_better > 5))

nrow(filter(selection_results,BT_dependent_Forest_better > 2))/nrow(selection_results)




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

drift_results$samplesize<-NA
samplesizes<-c(306,113,21,78,81,48,34,200,103,290,200,306,306,306,306,306,NA)

for (i in 1:length(unique(drift_results$Sample))) {
drift_results[drift_results$Sample %in% unique(drift_results$Sample)[i],]$samplesize<-samplesizes[i]
}

selection_results$samplesize<-NA

for (i in 1:length(unique(selection_results$Sample))) {
  selection_results[selection_results$Sample %in% unique(selection_results$Sample)[i],]$samplesize<-samplesizes[i]
}
  
#Horizontal transmission listing
drift_results$horizontaltransmission<-"No"
drift_results[grep("horizontal",drift_results$Sample),]$horizontaltransmission<-"Yes"

selection_results$horizontaltransmission<-"No"
selection_results[grep("horizontal",selection_results$Sample),]$horizontaltransmission<-"Yes"

#amount of diversity observed - something wrong here...
drift_results[drift_results$NumberOfVariantsInModel %in% 1,]$NumberOfVariantsInModel<-4
drift_results$diversityrepresented<-drift_results$NumberOfVariantsObserved/drift_results$NumberOfVariantsInModel

drift_results[drift_results$diversityrepresented %in% 1.5,]$NumberOfVariantsInModel<-4

selection_results$diversityrepresented<-selection_results$NumberOfVariantsObserved/selection_results$NumberOfVariantsInModel

drift_results<-drift_results[is.na(drift_results$Sample)==F,]
selection_results<-selection_results[is.na(selection_results$Sample)==F,]

drift_results$selectionwronglysupportedsignificant<-as.integer(drift_results$chisq_p_selectionwronglysupported<0.05)
drift_results$selectiontoonewronglysupportedsignificant<-as.integer(drift_results$chisq_p_selectiontoonestatewronglysupported<0.05)

levels(selection_results$Tree)<-c(levels(selection_results$Tree),"Grafentree","Earlytree","Grafentree","Latetree","Onetree")
selection_results[selection_results$Tree %in% "EarlyTree",]$Tree<-"Earlytree"
selection_results[selection_results$Tree %in% "Grafe",]$Tree<-"Grafentree"
selection_results[selection_results$Tree %in% "GrafenTree",]$Tree<-"Grafentree"
selection_results[selection_results$Tree %in% "LateTree",]$Tree<-"Latetree"
selection_results[selection_results$Tree %in% "OneTree",]$Tree<-"Onetree"



drift_results<-mutate(drift_results,BT_MultistateFree_better = 2*(drift_results[,]$BayestraitsMultistate_Free-drift_results[,]$BayestraitsMultistate_Equal) )
selection_results<-mutate(selection_results,BT_MultistateFree_better = 2*(selection_results[,]$BayestraitsMultistate_Free-selection_results[,]$BayestraitsMultistate_Equal) )

drift_results<-mutate(drift_results,chisq_BT_MultitstateFree_better=(1-pchisq(drift_results[,]$BT_MultistateFree_better, 10 )))
selection_results<-mutate(selection_results,chisq_BT_MultitstateFree_better=(1-pchisq(selection_results[,]$BT_MultistateFree_better, 10 )))

drift_results<-mutate(drift_results,chisq_BT_CladeAdependent_better=(1-pchisq(drift_results[,]$BT_dependent_CladeA_better, 4 )))
drift_results<-mutate(drift_results,chisq_BT_Forestdependent_better=(1-pchisq(drift_results[,]$BT_dependent_Forest_better, 4 )))

selection_results<-mutate(selection_results,chisq_BT_CladeAdependent_better=(1-pchisq(selection_results[,]$BT_dependent_CladeA_better, 4 )))
selection_results<-mutate(selection_results,chisq_BT_Forestdependent_better=(1-pchisq(selection_results[,]$BT_dependent_Forest_better, 4 )))


selection_results[selection_results$BayestraitsMultistate_Free %in% min(selection_results$BayestraitsMultistate_Free,na.rm=T),]$BayestraitsMultistate_Free<-NA
selection_results[selection_results$BayestraitsMultistate_Equal %in% min(selection_results$BayestraitsMultistate_Equal,na.rm=T),]$BayestraitsMultistate_Equal<-NA

drift_results[drift_results$BayestraitsMultistate_Free %in% min(drift_results$BayestraitsMultistate_Free,na.rm=T),]$BayestraitsMultistate_Free<-NA
drift_results[drift_results$BayestraitsMultistate_Equal %in% min(drift_results$BayestraitsMultistate_Equal,na.rm=T),]$BayestraitsMultistate_Equal<-NA


drift_results$BT_MultitstateFree_wronglysupported<-as.integer(drift_results$chisq_BT_MultitstateFree_better<0.05)
drift_results$BT_CladeAdependent_wronglysupported<-as.integer(drift_results$chisq_BT_CladeAdependent_better<0.05)
drift_results$BT_Forestdependent_wronglysupported<-as.integer(drift_results$chisq_BT_Forestdependent_better<0.05)

drift_results$X.1<-c(1:nrow(drift_results))

drift_results[filter(drift_results,(BayestraitsDependent_CladeA-BayestraitsIndependent_CladeA)>50)$X.1,]$BayestraitsDependent_CladeA<-NA

drift_results[filter(drift_results,(BayestraitsIndependent_CladeA-BayestraitsDependent_CladeA)>50)$X.1,]$BayestraitsIndependent_CladeA<-NA

drift_results[filter(drift_results,(BayestraitsIndependent_Ecology-BayestraitsDependent_Ecology)>5)$X.1,]$BayestraitsIndependent_Ecology<-NA

drift_results[filter(drift_results,(BayestraitsDependent_Ecology-BayestraitsIndependent_Ecology)>15 & BayestraitsIndependent_Ecology>"-220")$X.1,]$BayestraitsIndependent_Ecology<-NA


drift_results<-mutate(drift_results,BT_Multistate_nonconvergence=(drift_results[,]$BayestraitsMultistate_Free - drift_results[,]$BayestraitsMultistate_Equal<0))



selection_results$X.1<-c(1:nrow(selection_results))
selection_results[filter(selection_results,(BayestraitsDependent_CladeA-BayestraitsIndependent_CladeA)>50)$X.1,]$BayestraitsDependent_CladeA<-NA

selection_results[filter(selection_results,(BayestraitsIndependent_CladeA-BayestraitsDependent_CladeA)>50)$X.1,]$BayestraitsIndependent_CladeA<-NA

selection_results[filter(selection_results,(BayestraitsIndependent_Forest-BayestraitsDependent_Forest)>5)$X.1,]$BayestraitsIndependent_Forest<-NA
selection_results$X.1<-c(1:nrow(selection_results))
selection_results[filter(selection_results,(BayestraitsDependent_Forest-BayestraitsIndependent_Forest)>20)$X.1,]$BayestraitsIndependent_Forest<-NA

selection_results<-mutate(selection_results,BT_Multistate_nonconvergence=(selection_results[,]$BayestraitsMultistate_Free - selection_results[,]$BayestraitsMultistate_Equal<0))

selection_results$BT_MultitstateFree_correctlysupported<-as.integer(selection_results$chisq_BT_MultitstateFree_better<0.05)
selection_results$BT_CladeAdependent_wronglysupported<-as.integer(selection_results$chisq_BT_CladeAdependent_better<0.05)
selection_results$BT_Forestdependent_wronglysupported<-as.integer(selection_results$chisq_BT_Forestdependent_better<0.05)



#write new files to csv
write.csv(drift_results,file="CulturalPhylogenetics_SimulationResults_Drift_ForAnalysis.csv")
write.csv(selection_results,file="CulturalPhylogenetics_SimulationResults_Selection_ForAnalysis.csv")








library(ape)
library(geiger)
library(phytools)
library(OUwie)
library(dplyr)
library(btw)

setwd("~/ownCloud/Documents/CulturalPhylogenetics/Code/PamaNyungan/")

drift_results<-read.csv("CulturalPhylogenetics_SimulationResults_Drift_ForAnalysis.csv")
selection_results<-read.csv("CulturalPhylogenetics_SimulationResults_Selection_ForAnalysis.csv")


[1] "X.2"                                         "X.1"                                        
[3] "X"                                           "Tree"                                       
[5] "NumberOfVariantsInModel"                     "RateOfChange"                               
[7] "Repetition"                                  "Sample"                                     
[9] "NumberOfVariantsObserved"                    "loglik_drift"                               
[11] "loglik_selectiononestate"                    "loglik_pathways"                            
[13] "loglik_pathwaystoone"                        "loglik_unconstrained"                       
[15] "phylosigLambda"                              "phylosigK"                                  
[17] "BayestraitsIndependent_CladeA"               "BayestraitsDependent_CladeA"                
[19] "BayestraitsIndependent_Forest"               "BayestraitsDependent_Forest"                
[21] "loglik_pathwaysfromone"                      "BayestraitsMultistate_Free"                 
[23] "BayestraitsMultistate_Equal"                 "Pathway_impossible"                         
[25] "simulation_variant"                          "ml_model_supported"                         
[27] "chisq_p_selectionwronglysupported"           "chisq_p_selectiontoonestatewronglysupported"
[29] "BT_dependent_CladeA_better"                  "BT_dependent_Forest_better"                 
[31] "samplesize"                                  "horizontaltransmission"                     
[33] "diversityrepresented"                        "selectionwronglysupportedsignificant"       
[35] "selectiontoonewronglysupportedsignificant"   "BT_MultistateFree_better"                   
[37] "chisq_BT_MultitstateFree_better"             "chisq_BT_CladeAdependent_better"            
[39] "chisq_BT_Forestdependent_better"  


[1] "X.2"                                           "X.1"                                          
[3] "X"                                             "Tree"                                         
[5] "NumberOfVariantsInModel"                       "RateOfChange"                                 
[7] "Repetition"                                    "Sample"                                       
[9] "NumberOfVariantsObserved"                      "loglik_drift"                                 
[11] "loglik_selectiononestate"                      "loglik_pathways"                              
[13] "loglik_pathwaystoone"                          "loglik_unconstrained"                         
[15] "phylosigLambda"                                "phylosigK"                                    
[17] "BayestraitsIndependent_CladeA"                 "BayestraitsDependent_CladeA"                  
[19] "BayestraitsIndependent_Forest"                 "BayestraitsDependent_Forest"                  
[21] "loglik_pathwaysfromone"                        "BayestraitsMultistate_Free"                   
[23] "BayestraitsMultistate_Equal"                   "Pathway_impossible"                           
[25] "simulation_variant"                            "ml_model_supported"                           
[27] "chisq_p_selectioncorrectlysupported"           "chisq_p_selectiontoonestatecorrectlysupported"
[29] "BT_dependent_CladeA_better"                    "BT_dependent_Forest_better"                   
[31] "samplesize"                                    "horizontaltransmission"                       
[33] "diversityrepresented"                          "selectioncorrectlysupportedsignificant"       
[35] "selectiontoonecorrectlysupportedsignificant"   "BT_MultistateFree_better"                     
[37] "chisq_BT_MultitstateFree_better"               "chisq_BT_CladeAdependent_better"              
[39] "chisq_BT_Forestdependent_better" 



#Analysis of the effect of sample size


samplesize_analysis<-matrix(nrow=length(unique(drift_results$Sample)),ncol=2)
samplesize_analysis<-as.data.frame(samplesize_analysis)
colnames(samplesize_analysis)<-c("Sample","samplesize")
samplesize_analysis$samplesize<-c(172,102,28,42,36,27,39,129,86,167,108,127,172,172,172,172)
samplesize_analysis$Sample <-unique(drift_results$Sample)

diversity_drift<-drift_results %>% group_by(Sample) %>% summarise(mean(diversityrepresented))
diversity_drift<-as.data.frame(diversity_drift)
samplesize_analysis<-full_join(samplesize_analysis,diversity_drift,by="Sample")
colnames(samplesize_analysis)<-c(colnames(samplesize_analysis)[1:2],"drift_diversityrepresented")

diversity_selection<-selection_results %>% group_by(Sample) %>% summarise(mean(diversityrepresented))
diversity_selection<-as.data.frame(diversity_selection)
samplesize_analysis<-full_join(samplesize_analysis,diversity_selection,by="Sample")
colnames(samplesize_analysis)<-c(colnames(samplesize_analysis)[1:3],"selection_diversityrepresented")

drift_iterations<-drift_results %>% group_by(Sample) %>% summarise(n_distinct(X))
drift_iterations<-as.data.frame(drift_iterations)
samplesize_analysis<-full_join(samplesize_analysis,drift_iterations,by="Sample")
colnames(samplesize_analysis)<-c(colnames(samplesize_analysis)[1:4],"drift_iterations")

selection_iterations<-selection_results %>% group_by(Sample) %>% summarise(n_distinct(X))
selection_iterations<-as.data.frame(selection_iterations)
samplesize_analysis<-full_join(samplesize_analysis,selection_iterations,by="Sample")
colnames(samplesize_analysis)<-c(colnames(samplesize_analysis)[1:5],"selection_iterations")

drift_selectionwronglysupported<-drift_results %>% filter(chisq_p_selectionwronglysupported_df<0.05) %>% group_by(Sample) %>% summarise(n_distinct(X))
drift_selectionwronglysupported<-as.data.frame(drift_selectionwronglysupported)
samplesize_analysis<-full_join(samplesize_analysis,drift_selectionwronglysupported,by="Sample")
colnames(samplesize_analysis)<-c(colnames(samplesize_analysis)[1:6],"drift_selectionwronglysupported")

samplesize_analysis$proportion_drift_selectionwronglysupported<-samplesize_analysis$drift_selectionwronglysupported/samplesize_analysis$drift_iterations
colnames(samplesize_analysis)<-c(colnames(samplesize_analysis)[1:7],"proportion_drift_selectionwronglysupported")

selection_selectionwronglynotsupported<-selection_results %>% filter(chisq_p_selectionwcorrectlysupported>0.05) %>% group_by(Sample) %>% summarise(n_distinct(X))
selection_selectionwronglynotsupported<-as.data.frame(selection_selectionwronglynotsupported)
samplesize_analysis<-full_join(samplesize_analysis,selection_selectionwronglynotsupported,by="Sample")
colnames(samplesize_analysis)<-c(colnames(samplesize_analysis)[1:8],"selection_selectionwronglynotsupported")

samplesize_analysis$proportion_selection_selectionwronglynotsupported<-samplesize_analysis$selection_selectionwronglynotsupported/samplesize_analysis$selection_iterations
colnames(samplesize_analysis)<-c(colnames(samplesize_analysis)[1:9],"proportion_selection_selectionwronglynotsupported")

selection_selectiontoonewronglynotsupported<-selection_results %>% filter(chisq_p_selectiontoonestatecorrectlysupported>0.05) %>% group_by(Sample) %>% summarise(n_distinct(X))
selection_selectiontoonewronglynotsupported<-as.data.frame(selection_selectiontoonewronglynotsupported)
samplesize_analysis<-full_join(samplesize_analysis,selection_selectiontoonewronglynotsupported,by="Sample")
colnames(samplesize_analysis)<-c(colnames(samplesize_analysis)[1:10],"selection_selectiontoonewronglynotsupported")

samplesize_analysis$proportion_selection_selectiontoonewronglynotsupported<-samplesize_analysis$selection_selectiontoonewronglynotsupported/samplesize_analysis$selection_iterations
colnames(samplesize_analysis)<-c(colnames(samplesize_analysis)[1:11],"proportion_selection_selectiontoonewronglynotsupported")










drift_results %>% filter(horizontaltransmission %in% "No") %>% summarise(n_distinct(X))

drift_results %>% filter(chisq_p_selectionwronglysupported<0.05,horizontaltransmission %in% "No") %>% summarise(n_distinct(X))

drift_results %>% filter(chisq_p_selectionwronglysupported<0.05,horizontaltransmission %in% "No") %>% group_by(RateOfChange) %>% summarise(n_distinct(X))
drift_results %>% filter(horizontaltransmission %in% "No") %>% group_by(RateOfChange) %>% summarise(n_distinct(X))


drift_results %>% filter(chisq_p_selectionwronglysupported<0.05,horizontaltransmission %in% "No") %>% group_by(NumberOfVariantsObserved) %>% summarise(n_distinct(X))
drift_results %>% filter(horizontaltransmission %in% "No") %>% group_by(NumberOfVariantsObserved) %>% summarise(n_distinct(X))


drift_results %>% filter(chisq_p_selectionwronglysupported<0.05,horizontaltransmission %in% "Yes") %>% group_by(NumberOfVariantsObserved) %>% summarise(n_distinct(X))
drift_results %>% filter(horizontaltransmission %in% "Yes") %>% group_by(NumberOfVariantsObserved) %>% summarise(n_distinct(X))

drift_results %>% filter(chisq_p_selectionwronglysupported<0.05,horizontaltransmission %in% "Yes") %>% group_by(RateOfChange) %>% summarise(n_distinct(X))
drift_results %>% filter(horizontaltransmission %in% "Yes") %>% group_by(RateOfChange) %>% summarise(n_distinct(X))


drift_results %>% filter(chisq_p_selectionwronglysupported<0.05,horizontaltransmission %in% "No") %>% group_by(Tree) %>% summarise(n_distinct(X))
drift_results %>% filter(horizontaltransmission %in% "No") %>% group_by(Tree) %>% summarise(n_distinct(X))


drift_results %>% filter(chisq_p_selectionwronglysupported<0.05,horizontaltransmission %in% "Yes") %>% group_by(Tree) %>% summarise(n_distinct(X))
drift_results %>% filter(horizontaltransmission %in% "Yes") %>% group_by(Tree) %>% summarise(n_distinct(X))



nonconverged<-drift_results %>%  filter(BT_Multistate_nonconvergence %in% "TRUE") %>% group_by(Sample) %>% summarise(n_distinct(X))
converged<-drift_results %>%  filter(BT_Multistate_nonconvergence %in% "FALSE") %>% group_by(Sample) %>% summarise(n_distinct(X))


drift_results %>% filter(chisq_BT_MultitstateFree_better>0.05,horizontaltransmission %in% "Yes") %>% group_by(Sample) %>% summarise(n_distinct(X))
drift_results %>% filter(chisq_BT_MultitstateFree_better<0.05,horizontaltransmission %in% "Yes") %>% group_by(Sample) %>% summarise(n_distinct(X))

drift_results %>% filter(chisq_BT_MultitstateFree_better>0.05,horizontaltransmission %in% "No") %>% group_by(NumberOfVariantsObserved) %>% summarise(n_distinct(X))
drift_results %>% filter(chisq_BT_MultitstateFree_better<0.05,horizontaltransmission %in% "No") %>% group_by(NumberOfVariantsObserved) %>% summarise(n_distinct(X))


drift_results %>% filter(chisq_BT_MultitstateFree_better>0.05,horizontaltransmission %in% "Yes",BT_Multistate_nonconvergence %in% "FALSE") %>% group_by(Sample) %>% summarise(n_distinct(X))
drift_results %>% filter(chisq_BT_MultitstateFree_better<0.05,horizontaltransmission %in% "Yes",BT_Multistate_nonconvergence %in% "FALSE") %>% group_by(Sample) %>% summarise(n_distinct(X))


drift_results %>% filter(chisq_BT_MultitstateFree_better>0.05,horizontaltransmission %in% "No",BT_Multistate_nonconvergence %in% "FALSE") %>% group_by(NumberOfVariantsObserved) %>% summarise(n_distinct(X))
drift_results %>% filter(chisq_BT_MultitstateFree_better<0.05,horizontaltransmission %in% "No",BT_Multistate_nonconvergence %in% "FALSE") %>% group_by(NumberOfVariantsObserved) %>% summarise(n_distinct(X))

drift_results %>% filter(chisq_BT_MultitstateFree_better>0.05, BT_Multistate_nonconvergence %in% "FALSE") %>% group_by(Tree) %>% summarise(n_distinct(X))
drift_results %>% filter(chisq_BT_MultitstateFree_better<0.05, BT_Multistate_nonconvergence %in% "FALSE") %>% group_by(Tree) %>% summarise(n_distinct(X))


drift_results %>% filter(chisq_BT_CladeAdependent_better>0.05) %>% group_by(Tree) %>% summarise(n_distinct(X))
drift_results %>% filter(chisq_BT_CladeAdependent_better<0.05) %>% group_by(Tree) %>% summarise(n_distinct(X))

drift_results %>% filter(chisq_BT_CladeAdependent_better<0.05) %>% group_by(chisq_BT_CladeAdependent_better<0.05) %>% summarise(mean(phylosigK))


drift_results %>% filter(chisq_BT_Forestdependent_better>0.05) %>% group_by(Tree) %>% summarise(n_distinct(X))
drift_results %>% filter(chisq_BT_Forestdependent_better<0.05) %>% group_by(Tree) %>% summarise(n_distinct(X))


drift_results %>% filter(chisq_BT_Forestdependent_better<0.05) %>% group_by(Sample) %>% summarise(n_distinct(X))
drift_results %>% filter(chisq_BT_Forestdependent_better>0.05) %>% group_by(Sample) %>% summarise(n_distinct(X))


drift_results %>%  group_by(chisq_BT_MultitstateFree_better<0.5) %>% summarise(mean(phylosigK))




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

