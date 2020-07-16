#------------------------------------------------------------------------------------------
# Simulations for the manuscript - drift
#------------------------------------------------------------------------------------------

#Load necessary libraries

library(ape)
library(geiger)
library(phytools)
library(OUwie)
library(dplyr)
library(btw)


#------------------------------------------------------------------------------------------
# Load relevant WNAI data from GitHub
WNAIdata<-read_csv(url("https://raw.githubusercontent.com/dieterlukas/CulturalMacroevolution_Simulation/master/WNAI_data_forsimulation.csv"))
WNAIdata<-data.frame(WNAIdata)
data<-WNAIdata


#------------------------------------------------------------------------------------------

#Load tree constructed from the language classifications
Americantree<-read.nexus(url("https://raw.githubusercontent.com/dieterlukas/CulturalMacroevolution_Simulation/master/WNAI_tree_forsimulation.nex"))


#------------------------------------------------------------------------------------------
#Modify the tree for the analyses

#Add branch lengths to the tree based on Grafen's method (all tips equidistant from root, branch length depends onnumber of nodes between root and tip)
Grafentree<-compute.brlen(Americantree,method="Grafen")

#Add branch lengths to the tree assuming that all branches have the same length of one
Onetree<-compute.brlen(Americantree,1)

#Add branch lengths to the tree with early diversification and long branches to the tips
Earlytree<-compute.brlen(Americantree,method="Grafen",power=0.25)

#Add branch lengths to the tree with lots of recent diversification and long branches between clades
Latetree<-compute.brlen(Americantree,method="Grafen",power=1.5)

#Some analyses need a rooted, fully bifurcating tree
Grafentree<-root(Grafentree,node=173)
Grafentree<-multi2di(Grafentree)
Grafentree<-compute.brlen(Grafentree,method="Grafen")

Onetree<-root(Onetree,node=173)
Onetree<-multi2di(Onetree)
Onetree<-compute.brlen(Onetree,1)

Latetree<-root(Latetree,node=173)
Latetree<-multi2di(Latetree)
Latetree<-compute.brlen(Latetree,method="Grafen",power=1.5)

Earlytree<-root(Earlytree,node=173)
Earlytree<-multi2di(Earlytree)
Earlytree<-compute.brlen(Earlytree,method="Grafen",power=0.25)

#------------------------------------------------------------------------------------------
#Identify the major subclades in the tree
mrca(Americantree)["92","110"]  #192
mrca(Americantree)["3","1"] #174
mrca(Americantree)["172","72"] #240
mrca(Americantree)["6","160"] #193
mrca(Americantree)["44","119"] #227
mrca(Americantree)["162","10"] #211

#Create subsets of trees according to the major language groups
cladeA <- tips(Americantree, 192) # 102 societies "NorthernAmerind"
cladeB <- tips(Americantree, 174) # 28 societies "NaDene"
cladeC <- tips(Americantree, 240) # 42 societies "CentralAmerind"
cladeD <- tips(Americantree, 193) # 36 societies "Penutian"
cladeE <- tips(Americantree, 227) # 27 societies "Hokan"
cladeF <- tips(Americantree, 211) # 39 societies "Almosan"

cladeAtree<-keep.tip(Americantree,cladeA)
cladeBtree<-keep.tip(Americantree,cladeB)
cladeCtree<-keep.tip(Americantree,cladeC)
cladeDtree<-keep.tip(Americantree,cladeD)
cladeEtree<-keep.tip(Americantree,cladeE)
cladeFtree<-keep.tip(Americantree,cladeF)


#For the analyses asking whether selection operated differently in different clades, 
#we split the tree according to the largest clades - CladeA versus the remaining societies

#Create dummy variables designating subclade membership for each society 
data$cladeA<-rep(0,172)
data[data$tribes %in% cladeA,]$cladeA<-1
cladeAmember<-data$cladeA
names(cladeAmember)<-data$tribes

#Label all the nodes of the respective subclades
simmapattemptCladeA<-make.simmap(Grafentree,cladeAmember)

#We also want that selection operates differently in one of the ecoregions - all of them that are forests
data$Forests<-as.numeric(grepl("Forest",data$ecoregion))
Forestmember<-data$Forests
names(Forestmember)<-data$tribes

#Label all the nodes according to whether ancestral society was likely to be associated with forests
simmapattemptForests<-make.simmap(Grafentree,Forestmember)



#------------------------------------------------------------------------------------------
#Create the matrices for the simulations of the data on the phylogenetic tree

#Matrices reflecting drift, where transitions between all variants are equally likely

# Two states (labeled a,b); slow rate
drift_two_slow<-matrix(c(-0.2,0.2,  0.2,-0.2),2,2,dimnames=list(letters[1:2],letters[1:2]))

# Two states (a,b); medium rate
drift_two_medium<-matrix(c(-0.8,0.8,  0.8,-0.8),2,2,dimnames=list(letters[1:2],letters[1:2]))

# Two states (a,b); fast rate
drift_two_fast<-matrix(c(-1.4,1.4,  1.4,-1.4),2,2,dimnames=list(letters[1:2],letters[1:2]))


# Four states (a,b,c,d); slow rate
drift_four_slow<-matrix(c(-0.2,0.2,0.2,0.2, 0.2,-0.2,0.2,0.2, 0.2,0.2,-0.2,0.2, 0.2,0.2,0.2,-0.2),4,4,dimnames=list(letters[1:4],letters[1:4]))

# Four states (a,b,c,d); medium rate
drift_four_medium<-matrix(c(-0.8,0.8,0.8,0.8, 0.8,-0.8,0.8,0.8, 0.8,0.8,-0.8,0.8, 0.8,0.8,0.8,-0.8),4,4,dimnames=list(letters[1:4],letters[1:4]))

# Four states (a,b,c,d); fast rate
drift_four_fast<-matrix(c(-1.4,1.4,1.4,1.4, 1.4,-1.4,1.4,1.4, 1.4,1.4,-1.4,1.4, 1.4,1.4,1.4,-1.4),4,4,dimnames=list(letters[1:4],letters[1:4]))




#------------------------------------------------------------------------------------------
#Set up the matrices for the phylogenetic reconstructions

#In some simulations and especially in some subsets, we will not find all four variants so we need to make this option available
#The numbers assign which transition rates are supposed to be zero, the same as other transitions, or unique

#We first assume that transitions to and from one state are differently likely than other transitions
reconstruction_selection_four <- matrix (c(0,1,1,1, 2,0,3,3, 2,3,0,3, 2,3,3,0), nrow=4)
reconstruction_selection_three <- matrix (c(0,1,1, 2,0,3, 2,3,0), nrow=3)
reconstruction_selection_two <- matrix (c(0,1, 2,0), nrow=2)
reconstruction_selection_one <- matrix (c(1), nrow=1)

#Next we assume that transitions have to occur along paths from a <-> b <-> c <-> d
reconstruction_selection_four_pathway <- matrix (c(0,1,0,0, 2,0,3,0, 0,4,0,5, 0,0,6,0), nrow=4)
reconstruction_selection_three_pathway <- matrix (c(0,1,0, 2,0,3, 0,4,0), nrow=3)
reconstruction_selection_two_pathway <- matrix (c(0,1, 2,0), nrow=2)
reconstruction_selection_one_pathway <- matrix (c(1), nrow=1)

#We assume that transitions have to occur along paths from a <-> b <-> c <-> d, with transitions to a occurring at a different rate
reconstruction_selection_four_pathway_toone <- matrix (c(0,1,1,1, 0,0,2,2, 0,2,0,2, 0,2,2,0), nrow=4)
reconstruction_selection_three_pathway_toone <- matrix (c(0,1,1, 0,0,2, 0,2,0), nrow=3)
reconstruction_selection_two_pathway_toone <- matrix (c(0,1, 0,0), nrow=2)
reconstruction_selection_one_pathway_toone <- matrix (c(1), nrow=1)

#We assume that transitions have to occur along paths from a <-> b <-> c <-> d, with transitions from a occurring at a different rate
reconstruction_selection_four_pathway_fromone <- matrix (c(0,0,0,0, 1,0,2,2, 1,2,0,2, 1,2,2,0), nrow=4)
reconstruction_selection_three_pathway_fromone <- matrix (c(0,0,0, 1,0,2, 1,2,0), nrow=3)
reconstruction_selection_two_pathway_fromone <- matrix (c(0,0, 1,0), nrow=2)
reconstruction_selection_one_pathway_fromone <- matrix (c(1), nrow=1)



#-------------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------


#The following sets up loops to repeat simulations across all four trees and using the different drift models

#Run the simulation generating tip data under a drift model, testing whether the reconstruction wrongly favours that the data fit a model of selection

tree_variants <- c("Grafentree","Onetree","Earlytree","Latetree")

drift_variants<-c("drift_two_slow","drift_two_medium","drift_two_fast","drift_four_slow","drift_four_medium","drift_four_fast")

counter<-1

oldw <- getOption("warn")
options(warn = -1)

#Change the number of repetitions to modify the number of indepedent simulations with a given drift model on a given tree are being performed and analysed
repetitions<-50

#We prepare a data frame to store all the output
#Each row will have the results from a single simulation (with a given simulation model on a given tree)
#The columns contain basic information on the settings and the output from the various reconstruction methods
drift_results<-matrix(data=NA,nrow=500000,ncol=21)
colnames(drift_results)<-c("Tree","NumberOfVariantsInModel","RateOfChange","Repetition","Sample","NumberOfVariantsObserved","loglik_drift","loglik_selectiononestate","loglik_pathways","loglik_pathwaystoone","loglik_unconstrained","phylosigLambda","phylosigK","BayestraitsIndependent_CladeA","BayestraitsDependent_CladeA","BayestraitsIndependent_Forest","BayestraitsDependent_Forest","loglik_pathwaysfromone","BayestraitsMultistate_Free","BayestraitsMultistate_Equal","Pathway_impossible")


#This will start the loop; 
#because of the large number of analyses for each simulation (different techniques, various subsamples), this takes noticeable computer time

for (tree_variant in 1:length(tree_variants)) {
  
  tree_used<-tree_variants[tree_variant]
  if(tree_used=="Grafentree") currenttree<-Grafentree
  if(tree_used=="Onetree") currenttree<-Onetree
  if(tree_used=="Earlytree") currenttree<-Earlytree
  if(tree_used=="Latetree") currenttree<-Latetree
  
  for (drift_variant in 1:length(drift_variants)) {
    
    driftvariant_used<-drift_variants[drift_variant]
    if(driftvariant_used=="drift_two_slow") current_driftvariant<-drift_two_slow
    if(driftvariant_used=="drift_two_medium") current_driftvariant<-drift_two_medium
    if(driftvariant_used=="drift_two_fast") current_driftvariant<-drift_two_fast
    if(driftvariant_used=="drift_four_slow") current_driftvariant<-drift_four_slow
    if(driftvariant_used=="drift_four_medium") current_driftvariant<-drift_four_medium
    if(driftvariant_used=="drift_four_fast") current_driftvariant<-drift_four_fast
    
    for (repetition in 1:repetitions) {
      
      
      #Simulate the data with the specified tree and the specified drift model
      
      simulatedtips_neutral<-sim.Mk(tree=currenttree,Q=current_driftvariant)
      
      #store the data in case something breaks
      simulatedtips_problem<-simulatedtips_neutral
      
      
      #Convert the data for the analyses
      #The matrix for the selection reconstruction assumes that the first value is the one under selection
      #We are making the assumption that if there is selection, it would have acted on the variant present among the largest number of societies
      #We therefore assign the value that appears most commonly among the simulated tip values as the first value
      levels(simulatedtips_neutral)<-c(levels(simulatedtips_neutral),"1","2","3","4")
      simulatedtips_neutral[simulatedtips_neutral==names(sort(table(simulatedtips_problem),decreasing=TRUE)[1])]<-1
      simulatedtips_neutral[simulatedtips_neutral==names(sort(table(simulatedtips_problem),decreasing=TRUE)[2])]<-2
      simulatedtips_neutral[simulatedtips_neutral==names(sort(table(simulatedtips_problem),decreasing=TRUE)[3])]<-3
      simulatedtips_neutral[simulatedtips_neutral==names(sort(table(simulatedtips_problem),decreasing=TRUE)[4])]<-4
      
      simulatedtips_neutral<-droplevels(simulatedtips_neutral,exclude=levels(simulatedtips_problem))
      
      if (length(table(simulatedtips_neutral))!=1){
        
        #Start filling in the data in the respective column
        drift_results[counter,1]<-tree_used
        ifelse(drift_variant<4,drift_results[counter,2]<-2,drift_results[counter,2]<-4)
        ifelse(drift_variant<4,drift_results[counter,3]<-drift_variant,drift_results[counter,3]<-drift_variant-3)
        drift_results[counter,4]<-repetition
        
        
        
        #Start the phylogenetic reconstruction with the simulated data
        
        
        #Start with the full sample of 172 societies
        
        drift_results[counter,5]<-"full"
        drift_results[counter,6]<-length(unique(simulatedtips_neutral))
        
        
        # the phylogenetic reconstruction with the drift model, the model under which the data were simulated
        equal_reconstruction <- ace(simulatedtips_neutral, currenttree, type="discrete", model="ER")
        
        drift_results[counter,7]<-equal_reconstruction$loglik
        
        
        
        # the phylogenetic reconstruction with a one-state selection model
        transitionmatrix<-reconstruction_selection_four
        if(length(table(simulatedtips_neutral))==3) transitionmatrix<-reconstruction_selection_three
        if(length(table(simulatedtips_neutral))==2) transitionmatrix<-reconstruction_selection_two
        if(length(table(simulatedtips_neutral))==1) transitionmatrix<-reconstruction_selection_one
        
        selection_to_a_reconstruction_twofold <- ace(simulatedtips_neutral, currenttree, type="discrete", model=transitionmatrix)
        
        drift_results[counter,8]<-selection_to_a_reconstruction_twofold$loglik
        
        
        
        # the phylogenetic reconstruction with a pathway selection model
        transitionmatrix<-reconstruction_selection_four_pathway
        if(length(table(simulatedtips_neutral))==3) transitionmatrix<-reconstruction_selection_three_pathway
        if(length(table(simulatedtips_neutral))==2) transitionmatrix<-reconstruction_selection_two_pathway
        if(length(table(simulatedtips_neutral))==1) transitionmatrix<-reconstruction_selection_one_pathway
        
        selection_to_a_reconstruction_pathway_slow <- fitDiscrete(dat=simulatedtips_neutral, phy=currenttree, model="meristic")
        
        drift_results[counter,9]<-selection_to_a_reconstruction_pathway_slow$opt$lnL
        
        
        # the phylogenetic reconstruction with an alternative pathway selection model
        transitionmatrix<-reconstruction_selection_four_pathway_toone
        if(length(table(simulatedtips_neutral))==3) transitionmatrix<-reconstruction_selection_three_pathway_toone
        if(length(table(simulatedtips_neutral))==2) transitionmatrix<-reconstruction_selection_two_pathway_toone
        if(length(table(simulatedtips_neutral))==1) transitionmatrix<-reconstruction_selection_one_pathway_toone
        
        selection_to_a_reconstruction_pathway_toone <- ace(simulatedtips_neutral, currenttree, type="discrete", model=transitionmatrix)
        
        drift_results[counter,10]<-selection_to_a_reconstruction_pathway_toone$loglik
        
        
        # the phylogenetic reconstruction with an alternative pathway selection model
        transitionmatrix<-reconstruction_selection_four_pathway_fromone
        if(length(table(simulatedtips_neutral))==3) transitionmatrix<-reconstruction_selection_three_pathway_fromone
        if(length(table(simulatedtips_neutral))==2) transitionmatrix<-reconstruction_selection_two_pathway_fromone
        if(length(table(simulatedtips_neutral))==1) transitionmatrix<-reconstruction_selection_one_pathway_fromone
        
        selection_to_a_reconstruction_pathway_fromone <- ace(simulatedtips_neutral, currenttree, type="discrete", model=transitionmatrix)
        
        drift_results[counter,18]<-selection_to_a_reconstruction_pathway_fromone$loglik
        
        
        
        
        
        
        # the phylogenetic reconstruction with an unconstrained model
        unconstrained_reconstruction <- ace(simulatedtips_neutral, currenttree, type="discrete", model="ARD")
        
        drift_results[counter,11]<-unconstrained_reconstruction$loglik
        
        resultphylosiglambda<-phylosig(currenttree,as.numeric(simulatedtips_neutral),method="lambda") 
        drift_results[counter,12]<-resultphylosiglambda$lambda
        resultphylosigK<-phylosig(currenttree,as.numeric(simulatedtips_neutral),method="K") 
        drift_results[counter,13]<-resultphylosigK[1]
        
        
        #For the simulations where there are only two variants, we check whether by chance the variants
        #end up isolated in Clade A or associated with societies living in forest ecoregions
        #For this we use Bayestraits Discrete, assessing the likelihood of independent versus dependent evolution
        if(length(table(simulatedtips_neutral))==2) {
          
          bayestraitsdiscretedata<-matrix(NA,nrow=length(simulatedtips_neutral),ncol=3)
          colnames(bayestraitsdiscretedata)<-c("Species","Simulated","CladeA")
          rownames(bayestraitsdiscretedata)<-c(1:length(simulatedtips_neutral))
          bayestraitsdiscretedata[,1]<-names(simulatedtips_neutral)
          bayestraitsdiscretedata[,2]<-as.numeric(as.character(simulatedtips_neutral))
          bayestraitsdiscretedata[,2]<-as.numeric(bayestraitsdiscretedata[,2])-as.numeric(min(bayestraitsdiscretedata[,2]))
          bayestraitsdiscretedata[bayestraitsdiscretedata[,2]!=0,2]<-1
          bayestraitsdiscretedata[,3]<- as.numeric(as.character(cladeAmember))
          bayestraitsdiscretedata<-as.data.frame(bayestraitsdiscretedata)
          
          command_vec1 <- c("2", "1") #option 1 = 2 discrete independent; option 2 = 1 maximum likelihood
          results_1 <- bayestraits(bayestraitsdiscretedata, currenttree, command_vec1)
          log_1 <- results_1$Log
          drift_results[counter,14]<-log_1$results$Lh
          
          command_vec2 <- c("3", "1") #option 1 = 3 discrete dependent; option 2 = 1 maximum likelihood
          results_2 <- bayestraits(bayestraitsdiscretedata, currenttree, command_vec2)
          log_2 <- results_2$Log
          drift_results[counter,15]<-log_2$results$Lh
          
          
          bayestraitsdiscretedata<-matrix(NA,nrow=length(simulatedtips_neutral),ncol=3)
          colnames(bayestraitsdiscretedata)<-c("Species","Simulated","Forestmember")
          rownames(bayestraitsdiscretedata)<-c(1:length(simulatedtips_neutral))
          bayestraitsdiscretedata[,1]<-names(simulatedtips_neutral)
          bayestraitsdiscretedata[,2]<-as.numeric(as.character(simulatedtips_neutral))
          bayestraitsdiscretedata[,2]<-as.numeric(bayestraitsdiscretedata[,2])-as.numeric(min(bayestraitsdiscretedata[,2]))
          bayestraitsdiscretedata[bayestraitsdiscretedata[,2]!=0,2]<-1
          bayestraitsdiscretedata[,3]<- as.numeric(as.character(Forestmember))
          bayestraitsdiscretedata<-as.data.frame(bayestraitsdiscretedata)
          
          command_vec3 <- c("2", "1") #option 1 = 2 discrete independent; option 2 = 1 maximum likelihood
          results_3 <- bayestraits(bayestraitsdiscretedata, currenttree, command_vec3)
          log_3 <- results_3$Log
          drift_results[counter,16]<-log_3$results$Lh
          
          command_vec4 <- c("3", "1") #option 1 = 3 discrete dependent; option 2 = 1 maximum likelihood
          results_4 <- bayestraits(bayestraitsdiscretedata, currenttree, command_vec4)
          log_4 <- results_4$Log
          drift_results[counter,17]<-log_4$results$Lh
          
        } #end of bayestraits 2 option
        
        
        #For the simulations where there are more than two variants, we run the Bayestraits multistate option
        
        if(length(table(simulatedtips_neutral))>2) {
          
          bayestraitsdiscretedata<-matrix(NA,nrow=length(simulatedtips_neutral),ncol=2)
          colnames(bayestraitsdiscretedata)<-c("Species","Simulated")
          rownames(bayestraitsdiscretedata)<-c(1:length(simulatedtips_neutral))
          bayestraitsdiscretedata[,1]<-names(simulatedtips_neutral)
          bayestraitsdiscretedata[,2]<-as.numeric(as.character(simulatedtips_neutral))
          bayestraitsdiscretedata[,2]<-as.numeric(bayestraitsdiscretedata[,2])-1
          bayestraitsdiscretedata<-as.data.frame(bayestraitsdiscretedata)
          
          command_vec1 <- c("1", "1") #option 1 = 1 multistate; option 2 = 1 maximum likelihood
          results_1 <- bayestraits(bayestraitsdiscretedata, currenttree, command_vec1)
          log_1 <- results_1$Log
          drift_results[counter,19]<-log_1$results$Lh
          
          command_vec2 <- c("1", "1", "Restrict q01 q02 q03 q10 q12 q13 q20 q21 q23 q30 q31 q32") #option 1 = 1 multistate; option 2 = 1 maximum likelihood; option restrict all dependent
          results_2 <- bayestraits(bayestraitsdiscretedata, currenttree, command_vec2)
          log_2 <- results_2$Log
          drift_results[counter,20]<-log_2$results$Lh
          
          
        } #end of bayestraits multistate option
        
        
        
 
        
        
        counter<-counter+1
        
        
        
#--#--#--#--#--#--#        #--#--#--#--#--#--#
        
        #Repeat the analyses with the different subsamples
        
        
        # take the different subfamilies
        
        # Subset the tree and the data to match the respective clades
        
        cladeAsimulatedtips<-simulatedtips_neutral[names(simulatedtips_neutral) %in% cladeA]
        cladeAtree<-keep.tip(currenttree,cladeA)
        
        cladeAsimulatedtips<-droplevels(cladeAsimulatedtips)
        
        drift_results[counter,1]<-tree_used
        ifelse(drift_variant<4,drift_results[counter,2]<-2,drift_results[counter,2]<-4)
        ifelse(drift_variant<4,drift_results[counter,3]<-drift_variant,drift_results[counter,3]<-drift_variant-3)
        drift_results[counter,4]<-repetition
        
        drift_results[counter,5]<-"cladeA"
        drift_results[counter,6]<-length(unique(cladeAsimulatedtips))
        
        if (length(table(cladeAsimulatedtips))!=1) { 
          
          # the phylogenetic reconstruction with the drift model, the model under which the data were simulated
          equal_reconstruction <- ace(cladeAsimulatedtips, cladeAtree, type="discrete", model="ER")
          
          drift_results[counter,7]<-equal_reconstruction$loglik
          
          
          
          # the phylogenetic reconstruction with a slow one-state selection model
          transitionmatrix<-reconstruction_selection_four
          if(length(table(cladeAsimulatedtips))==3) transitionmatrix<-reconstruction_selection_three
          if(length(table(cladeAsimulatedtips))==2) transitionmatrix<-reconstruction_selection_two
          if(length(table(cladeAsimulatedtips))==1) transitionmatrix<-reconstruction_selection_one
          
          selection_to_a_reconstruction_twofold <- ace(cladeAsimulatedtips, cladeAtree, type="discrete", model=transitionmatrix)
          
          drift_results[counter,8]<-selection_to_a_reconstruction_twofold$loglik
          
          
          
          # the phylogenetic reconstruction with a pathway slow selection model
          transitionmatrix<-reconstruction_selection_four_pathway
          if(length(table(cladeAsimulatedtips))==3) transitionmatrix<-reconstruction_selection_three_pathway
          if(length(table(cladeAsimulatedtips))==2) transitionmatrix<-reconstruction_selection_two_pathway
          if(length(table(cladeAsimulatedtips))==1) transitionmatrix<-reconstruction_selection_one_pathway
          
          
          
          selection_to_a_reconstruction_pathway_slow <- fitDiscrete(dat=cladeAsimulatedtips, phy=cladeAtree, model="meristic")
          
          drift_results[counter,9]<-selection_to_a_reconstruction_pathway_slow$opt$lnL
          
          if(length(unique(cladeAsimulatedtips))>1) ifelse(names(table(cladeAsimulatedtips))[2]!=2, drift_results[counter,21]<-"TRUE",drift_results[counter,21]<-"FALSE")
          
          
          
          
          
          # the phylogenetic reconstruction with an alternative pathway selection model
          transitionmatrix<-reconstruction_selection_four_pathway_toone
          if(length(table(cladeAsimulatedtips))==3) transitionmatrix<-reconstruction_selection_three_pathway_toone
          if(length(table(cladeAsimulatedtips))==2) transitionmatrix<-reconstruction_selection_two_pathway_toone
          if(length(table(cladeAsimulatedtips))==1) transitionmatrix<-reconstruction_selection_one_pathway_toone
          
          selection_to_a_reconstruction_pathway_toone <- ace(cladeAsimulatedtips, cladeAtree, type="discrete", model=transitionmatrix)
          
          drift_results[counter,10]<-selection_to_a_reconstruction_pathway_toone$loglik
          
          # the phylogenetic reconstruction with an alternative pathway selection model
          transitionmatrix<-reconstruction_selection_four_pathway_fromone
          if(length(table(cladeAsimulatedtips))==3) transitionmatrix<-reconstruction_selection_three_pathway_fromone
          if(length(table(cladeAsimulatedtips))==2) transitionmatrix<-reconstruction_selection_two_pathway_fromone
          if(length(table(cladeAsimulatedtips))==1) transitionmatrix<-reconstruction_selection_one_pathway_fromone
          
          selection_to_a_reconstruction_pathway_fromone <- ace(cladeAsimulatedtips, cladeAtree, type="discrete", model=transitionmatrix)
          
          drift_results[counter,18]<-selection_to_a_reconstruction_pathway_fromone$loglik
          
          
          
          # the phylogenetic reconstruction with an unconstrained model
          unconstrained_reconstruction <- ace(cladeAsimulatedtips, cladeAtree, type="discrete", model="ARD")
          
          drift_results[counter,11]<-unconstrained_reconstruction$loglik
          
          resultphylosiglambda<-phylosig(cladeAtree,as.numeric(cladeAsimulatedtips),method="lambda") 
          drift_results[counter,12]<-resultphylosiglambda$lambda
          resultphylosigK<-phylosig(cladeAtree,as.numeric(cladeAsimulatedtips),method="K") 
          drift_results[counter,13]<-resultphylosigK[1]
          
        }
        
        counter<-counter+1
        
        
        
        
        #-#-#-#-#-#-#-#-#-#-#-#
        
        # Subset the tree and the data to match the respective clades
        
        cladeBsimulatedtips<-simulatedtips_neutral[names(simulatedtips_neutral) %in% cladeB]
        cladeBtree<-keep.tip(currenttree,cladeB)
        
        cladeBsimulatedtips<-droplevels(cladeBsimulatedtips)
        
        drift_results[counter,1]<-tree_used
        ifelse(drift_variant<4,drift_results[counter,2]<-2,drift_results[counter,2]<-4)
        ifelse(drift_variant<4,drift_results[counter,3]<-drift_variant,drift_results[counter,3]<-drift_variant-3)
        drift_results[counter,4]<-repetition
        
        drift_results[counter,5]<-"cladeB"
        drift_results[counter,6]<-length(unique(cladeBsimulatedtips))
        
        if (length(table(cladeBsimulatedtips))!=1) { 
          
          # the phylogenetic reconstruction with the drift model, the model under which the data were simulated
          equal_reconstruction <- ace(cladeBsimulatedtips, cladeBtree, type="discrete", model="ER")
          
          drift_results[counter,7]<-equal_reconstruction$loglik
          
          
          
          # the phylogenetic reconstruction with a slow one-state selection model
          transitionmatrix<-reconstruction_selection_four
          if(length(table(cladeBsimulatedtips))==3) transitionmatrix<-reconstruction_selection_three
          if(length(table(cladeBsimulatedtips))==2) transitionmatrix<-reconstruction_selection_two
          if(length(table(cladeBsimulatedtips))==1) transitionmatrix<-reconstruction_selection_one
          
          selection_to_a_reconstruction_twofold <- ace(cladeBsimulatedtips, cladeBtree, type="discrete", model=transitionmatrix)
          
          drift_results[counter,8]<-selection_to_a_reconstruction_twofold$loglik
          
          
          # the phylogenetic reconstruction with a pathway slow selection model
          transitionmatrix<-reconstruction_selection_four_pathway
          if(length(table(cladeBsimulatedtips))==3) transitionmatrix<-reconstruction_selection_three_pathway
          if(length(table(cladeBsimulatedtips))==2) transitionmatrix<-reconstruction_selection_two_pathway
          if(length(table(cladeBsimulatedtips))==1) transitionmatrix<-reconstruction_selection_one_pathway
          
          
          selection_to_a_reconstruction_pathway_slow <- fitDiscrete(dat=cladeBsimulatedtips, phy=cladeBtree, model="meristic")
          
          drift_results[counter,9]<-selection_to_a_reconstruction_pathway_slow$opt$lnL
          
          if(length(unique(cladeBsimulatedtips))>1) ifelse(names(table(cladeBsimulatedtips))[2]!=2, drift_results[counter,21]<-"TRUE",drift_results[counter,21]<-"FALSE")
          
          
          
          # the phylogenetic reconstruction with an alternative pathway selection model
          transitionmatrix<-reconstruction_selection_four_pathway_toone
          if(length(table(cladeBsimulatedtips))==3) transitionmatrix<-reconstruction_selection_three_pathway_toone
          if(length(table(cladeBsimulatedtips))==2) transitionmatrix<-reconstruction_selection_two_pathway_toone
          if(length(table(cladeBsimulatedtips))==1) transitionmatrix<-reconstruction_selection_one_pathway_toone
          
          selection_to_a_reconstruction_pathway_toone <- ace(cladeBsimulatedtips, cladeBtree, type="discrete", model=transitionmatrix)
          
          drift_results[counter,10]<-selection_to_a_reconstruction_pathway_toone$loglik
          
          # the phylogenetic reconstruction with an alternative pathway selection model
          transitionmatrix<-reconstruction_selection_four_pathway_fromone
          if(length(table(cladeBsimulatedtips))==3) transitionmatrix<-reconstruction_selection_three_pathway_fromone
          if(length(table(cladeBsimulatedtips))==2) transitionmatrix<-reconstruction_selection_two_pathway_fromone
          if(length(table(cladeBsimulatedtips))==1) transitionmatrix<-reconstruction_selection_one_pathway_fromone
          
          selection_to_a_reconstruction_pathway_fromone <- ace(cladeBsimulatedtips, cladeBtree, type="discrete", model=transitionmatrix)
          
          
          # the phylogenetic reconstruction with an unconstrained model
          unconstrained_reconstruction <- ace(cladeBsimulatedtips, cladeBtree, type="discrete", model="ARD")
          
          drift_results[counter,11]<-unconstrained_reconstruction$loglik
          
          resultphylosiglambda<-phylosig(cladeBtree,as.numeric(cladeBsimulatedtips),method="lambda") 
          drift_results[counter,12]<-resultphylosiglambda$lambda
          resultphylosigK<-phylosig(cladeBtree,as.numeric(cladeBsimulatedtips),method="K") 
          drift_results[counter,13]<-resultphylosigK[1]
          
        }
        
        counter<-counter+1
        
        
        
        
        #-#-#-#-#-#-#-#-#-#-#-#
        
        # Subset the tree and the data to match the respective clades
        
        cladeCsimulatedtips<-simulatedtips_neutral[names(simulatedtips_neutral) %in% cladeC]
        cladeCtree<-keep.tip(currenttree,cladeC)
        
        cladeCsimulatedtips<-droplevels(cladeCsimulatedtips)
        
        drift_results[counter,1]<-tree_used
        ifelse(drift_variant<4,drift_results[counter,2]<-2,drift_results[counter,2]<-4)
        ifelse(drift_variant<4,drift_results[counter,3]<-drift_variant,drift_results[counter,3]<-drift_variant-3)
        drift_results[counter,4]<-repetition
        
        drift_results[counter,5]<-"cladeC"
        drift_results[counter,6]<-length(unique(cladeCsimulatedtips))
        
        if (length(table(cladeCsimulatedtips))!=1) { 
          
          # the phylogenetic reconstruction with the drift model, the model under which the data were simulated
          equal_reconstruction <- ace(cladeCsimulatedtips, cladeCtree, type="discrete", model="ER")
          
          drift_results[counter,7]<-equal_reconstruction$loglik
          
          
          
          # the phylogenetic reconstruction with a slow one-state selection model
          transitionmatrix<-reconstruction_selection_four
          if(length(table(cladeCsimulatedtips))==3) transitionmatrix<-reconstruction_selection_three
          if(length(table(cladeCsimulatedtips))==2) transitionmatrix<-reconstruction_selection_two
          if(length(table(cladeCsimulatedtips))==1) transitionmatrix<-reconstruction_selection_one
          
          selection_to_a_reconstruction_twofold <- ace(cladeCsimulatedtips, cladeCtree, type="discrete", model=transitionmatrix)
          
          drift_results[counter,8]<-selection_to_a_reconstruction_twofold$loglik
          
          
          
          # the phylogenetic reconstruction with a pathway slow selection model
          transitionmatrix<-reconstruction_selection_four_pathway
          if(length(table(cladeCsimulatedtips))==3) transitionmatrix<-reconstruction_selection_three_pathway
          if(length(table(cladeCsimulatedtips))==2) transitionmatrix<-reconstruction_selection_two_pathway
          if(length(table(cladeCsimulatedtips))==1) transitionmatrix<-reconstruction_selection_one_pathway
          
          
          selection_to_a_reconstruction_pathway_slow <- fitDiscrete(dat=cladeCsimulatedtips, phy=cladeCtree, model="meristic")
          
          drift_results[counter,9]<-selection_to_a_reconstruction_pathway_slow$opt$lnL
          
          
          if(length(unique(cladeCsimulatedtips))>1) ifelse(names(table(cladeCsimulatedtips))[2]!=2, drift_results[counter,21]<-"TRUE",drift_results[counter,21]<-"FALSE")
          
          
          
          # the phylogenetic reconstruction with an alternative pathway selection model
          transitionmatrix<-reconstruction_selection_four_pathway_toone
          if(length(table(cladeCsimulatedtips))==3) transitionmatrix<-reconstruction_selection_three_pathway_toone
          if(length(table(cladeCsimulatedtips))==2) transitionmatrix<-reconstruction_selection_two_pathway_toone
          if(length(table(cladeCsimulatedtips))==1) transitionmatrix<-reconstruction_selection_one_pathway_toone
          
          selection_to_a_reconstruction_pathway_toone <- ace(cladeCsimulatedtips, cladeCtree, type="discrete", model=transitionmatrix)
          
          drift_results[counter,10]<-selection_to_a_reconstruction_pathway_toone$loglik
          
          # the phylogenetic reconstruction with an alternative pathway selection model
          transitionmatrix<-reconstruction_selection_four_pathway_fromone
          if(length(table(cladeCsimulatedtips))==3) transitionmatrix<-reconstruction_selection_three_pathway_fromone
          if(length(table(cladeCsimulatedtips))==2) transitionmatrix<-reconstruction_selection_two_pathway_fromone
          if(length(table(cladeCsimulatedtips))==1) transitionmatrix<-reconstruction_selection_one_pathway_fromone
          
          selection_to_a_reconstruction_pathway_fromone <- ace(cladeCsimulatedtips, cladeCtree, type="discrete", model=transitionmatrix)
          
          drift_results[counter,18]<-selection_to_a_reconstruction_pathway_fromone$loglik
          
          
          
          # the phylogenetic reconstruction with an unconstrained model
          unconstrained_reconstruction <- ace(cladeCsimulatedtips, cladeCtree, type="discrete", model="ARD")
          
          drift_results[counter,11]<-unconstrained_reconstruction$loglik
          
          resultphylosiglambda<-phylosig(cladeCtree,as.numeric(cladeCsimulatedtips),method="lambda") 
          drift_results[counter,12]<-resultphylosiglambda$lambda
          resultphylosigK<-phylosig(cladeCtree,as.numeric(cladeCsimulatedtips),method="K") 
          drift_results[counter,13]<-resultphylosigK[1]
          
        }
        
        counter<-counter+1
        
        
        
        
        #-#-#-#-#-#-#-#-#-#-#-#
        
        # Subset the tree and the data to match the respective clades
        
        cladeDsimulatedtips<-simulatedtips_neutral[names(simulatedtips_neutral) %in% cladeD]
        cladeDtree<-keep.tip(currenttree,cladeD)
        
        cladeDsimulatedtips<-droplevels(cladeDsimulatedtips)
        
        drift_results[counter,1]<-tree_used
        ifelse(drift_variant<4,drift_results[counter,2]<-2,drift_results[counter,2]<-4)
        ifelse(drift_variant<4,drift_results[counter,3]<-drift_variant,drift_results[counter,3]<-drift_variant-3)
        drift_results[counter,4]<-repetition
        
        drift_results[counter,5]<-"cladeD"
        drift_results[counter,6]<-length(unique(cladeDsimulatedtips))
        
        if (length(table(cladeDsimulatedtips))!=1) { 
          
          # the phylogenetic reconstruction with the drift model, the model under which the data were simulated
          equal_reconstruction <- ace(cladeDsimulatedtips, cladeDtree, type="discrete", model="ER")
          
          drift_results[counter,7]<-equal_reconstruction$loglik
          
          
          
          # the phylogenetic reconstruction with a slow one-state selection model
          transitionmatrix<-reconstruction_selection_four
          if(length(table(cladeDsimulatedtips))==3) transitionmatrix<-reconstruction_selection_three
          if(length(table(cladeDsimulatedtips))==2) transitionmatrix<-reconstruction_selection_two
          if(length(table(cladeDsimulatedtips))==1) transitionmatrix<-reconstruction_selection_one
          
          selection_to_a_reconstruction_twofold <- ace(cladeDsimulatedtips, cladeDtree, type="discrete", model=transitionmatrix)
          
          drift_results[counter,8]<-selection_to_a_reconstruction_twofold$loglik
          
          
          
          # the phylogenetic reconstruction with a pathway slow selection model
          transitionmatrix<-reconstruction_selection_four_pathway
          if(length(table(cladeDsimulatedtips))==3) transitionmatrix<-reconstruction_selection_three_pathway
          if(length(table(cladeDsimulatedtips))==2) transitionmatrix<-reconstruction_selection_two_pathway
          if(length(table(cladeDsimulatedtips))==1) transitionmatrix<-reconstruction_selection_one_pathway
          
          selection_to_a_reconstruction_pathway_slow <- fitDiscrete(dat=cladeDsimulatedtips, phy=cladeDtree, model="meristic")
          
          drift_results[counter,9]<-selection_to_a_reconstruction_pathway_slow$opt$lnL
          
          
          if(length(unique(cladeDsimulatedtips))>1) ifelse(names(table(cladeDsimulatedtips))[2]!=2, drift_results[counter,21]<-"TRUE",drift_results[counter,21]<-"FALSE")
          
          
          
          
          # the phylogenetic reconstruction with an alternative pathway selection model
          transitionmatrix<-reconstruction_selection_four_pathway_toone
          if(length(table(cladeDsimulatedtips))==3) transitionmatrix<-reconstruction_selection_three_pathway_toone
          if(length(table(cladeDsimulatedtips))==2) transitionmatrix<-reconstruction_selection_two_pathway_toone
          if(length(table(cladeDsimulatedtips))==1) transitionmatrix<-reconstruction_selection_one_pathway_toone
          
          selection_to_a_reconstruction_pathway_toone <- ace(cladeDsimulatedtips, cladeDtree, type="discrete", model=transitionmatrix)
          
          drift_results[counter,10]<-selection_to_a_reconstruction_pathway_toone$loglik
          
          # the phylogenetic reconstruction with an alternative pathway selection model
          transitionmatrix<-reconstruction_selection_four_pathway_fromone
          if(length(table(cladeDsimulatedtips))==3) transitionmatrix<-reconstruction_selection_three_pathway_fromone
          if(length(table(cladeDsimulatedtips))==2) transitionmatrix<-reconstruction_selection_two_pathway_fromone
          if(length(table(cladeDsimulatedtips))==1) transitionmatrix<-reconstruction_selection_one_pathway_fromone
          
          selection_to_a_reconstruction_pathway_fromone <- ace(cladeDsimulatedtips, cladeDtree, type="discrete", model=transitionmatrix)
          
          drift_results[counter,18]<-selection_to_a_reconstruction_pathway_fromone$loglik
          
          
          
          # the phylogenetic reconstruction with an unconstrained model
          unconstrained_reconstruction <- ace(cladeDsimulatedtips, cladeDtree, type="discrete", model="ARD")
          
          drift_results[counter,11]<-unconstrained_reconstruction$loglik
          
          resultphylosiglambda<-phylosig(cladeDtree,as.numeric(cladeDsimulatedtips),method="lambda") 
          drift_results[counter,12]<-resultphylosiglambda$lambda
          resultphylosigK<-phylosig(cladeDtree,as.numeric(cladeDsimulatedtips),method="K") 
          drift_results[counter,13]<-resultphylosigK[1]
          
        }
        
        counter<-counter+1
        
        
        
        
        #-#-#-#-#-#-#-#-#-#-#-#
        
        # Subset the tree and the data to match the respective clades
        
        cladeEsimulatedtips<-simulatedtips_neutral[names(simulatedtips_neutral) %in% cladeE]
        cladeEtree<-keep.tip(currenttree,cladeE)
        
        cladeEsimulatedtips<-droplevels(cladeEsimulatedtips)
        
        drift_results[counter,1]<-tree_used
        ifelse(drift_variant<4,drift_results[counter,2]<-2,drift_results[counter,2]<-4)
        ifelse(drift_variant<4,drift_results[counter,3]<-drift_variant,drift_results[counter,3]<-drift_variant-3)
        drift_results[counter,4]<-repetition
        
        drift_results[counter,5]<-"cladeE"
        drift_results[counter,6]<-length(unique(cladeEsimulatedtips))
        
        if (length(table(cladeEsimulatedtips))!=1) { 
          
          # the phylogenetic reconstruction with the drift model, the model under which the data were simulated
          equal_reconstruction <- ace(cladeEsimulatedtips, cladeEtree, type="discrete", model="ER")
          
          drift_results[counter,7]<-equal_reconstruction$loglik
          
          
          
          # the phylogenetic reconstruction with a slow one-state selection model
          transitionmatrix<-reconstruction_selection_four
          if(length(table(cladeEsimulatedtips))==3) transitionmatrix<-reconstruction_selection_three
          if(length(table(cladeEsimulatedtips))==2) transitionmatrix<-reconstruction_selection_two
          if(length(table(cladeEsimulatedtips))==1) transitionmatrix<-reconstruction_selection_one
          
          selection_to_a_reconstruction_twofold <- ace(cladeEsimulatedtips, cladeEtree, type="discrete", model=transitionmatrix)
          
          drift_results[counter,8]<-selection_to_a_reconstruction_twofold$loglik
          
          # the phylogenetic reconstruction with a pathway slow selection model
          transitionmatrix<-reconstruction_selection_four_pathway
          if(length(table(cladeEsimulatedtips))==3) transitionmatrix<-reconstruction_selection_three_pathway
          if(length(table(cladeEsimulatedtips))==2) transitionmatrix<-reconstruction_selection_two_pathway
          if(length(table(cladeEsimulatedtips))==1) transitionmatrix<-reconstruction_selection_one_pathway
          
          selection_to_a_reconstruction_pathway_slow <- fitDiscrete(dat=cladeEsimulatedtips, phy=cladeEtree, model="meristic")
          
          drift_results[counter,9]<-selection_to_a_reconstruction_pathway_slow$opt$lnL
          
          if(length(unique(cladeEsimulatedtips))>1) ifelse(names(table(cladeEsimulatedtips))[2]!=2, drift_results[counter,21]<-"TRUE",drift_results[counter,21]<-"FALSE")
          
          
          
          # the phylogenetic reconstruction with an alternative pathway selection model
          transitionmatrix<-reconstruction_selection_four_pathway_toone
          if(length(table(cladeEsimulatedtips))==3) transitionmatrix<-reconstruction_selection_three_pathway_toone
          if(length(table(cladeEsimulatedtips))==2) transitionmatrix<-reconstruction_selection_two_pathway_toone
          if(length(table(cladeEsimulatedtips))==1) transitionmatrix<-reconstruction_selection_one_pathway_toone
          
          selection_to_a_reconstruction_pathway_toone <- ace(cladeEsimulatedtips, cladeEtree, type="discrete", model=transitionmatrix)
          
          drift_results[counter,10]<-selection_to_a_reconstruction_pathway_toone$loglik
          
          # the phylogenetic reconstruction with an alternative pathway selection model
          transitionmatrix<-reconstruction_selection_four_pathway_fromone
          if(length(table(cladeEsimulatedtips))==3) transitionmatrix<-reconstruction_selection_three_pathway_fromone
          if(length(table(cladeEsimulatedtips))==2) transitionmatrix<-reconstruction_selection_two_pathway_fromone
          if(length(table(cladeEsimulatedtips))==1) transitionmatrix<-reconstruction_selection_one_pathway_fromone
          
          selection_to_a_reconstruction_pathway_fromone <- ace(cladeEsimulatedtips, cladeEtree, type="discrete", model=transitionmatrix)
          
          drift_results[counter,18]<-selection_to_a_reconstruction_pathway_fromone$loglik
          
          
          
          # the phylogenetic reconstruction with an unconstrained model
          unconstrained_reconstruction <- ace(cladeEsimulatedtips, cladeEtree, type="discrete", model="ARD")
          
          drift_results[counter,11]<-unconstrained_reconstruction$loglik
          
          resultphylosiglambda<-phylosig(cladeEtree,as.numeric(cladeEsimulatedtips),method="lambda") 
          drift_results[counter,12]<-resultphylosiglambda$lambda
          resultphylosigK<-phylosig(cladeEtree,as.numeric(cladeEsimulatedtips),method="K") 
          drift_results[counter,13]<-resultphylosigK[1]
          
        }
        
        counter<-counter+1
        
        
        
        
        #-#-#-#-#-#-#-#-#-#-#-#
        
        # Subset the tree and the data to match the respective clades
        
        cladeFsimulatedtips<-simulatedtips_neutral[names(simulatedtips_neutral) %in% cladeF]
        cladeFtree<-keep.tip(currenttree,cladeF)
        
        cladeFsimulatedtips<-droplevels(cladeFsimulatedtips)
        
        drift_results[counter,1]<-tree_used
        ifelse(drift_variant<4,drift_results[counter,2]<-2,drift_results[counter,2]<-4)
        ifelse(drift_variant<4,drift_results[counter,3]<-drift_variant,drift_results[counter,3]<-drift_variant-3)
        drift_results[counter,4]<-repetition
        
        drift_results[counter,5]<-"cladeF"
        drift_results[counter,6]<-length(unique(cladeFsimulatedtips))
        
        if (length(table(cladeFsimulatedtips))!=1) { 
          
          # the phylogenetic reconstruction with the drift model, the model under which the data were simulated
          equal_reconstruction <- ace(cladeFsimulatedtips, cladeFtree, type="discrete", model="ER")
          
          drift_results[counter,7]<-equal_reconstruction$loglik
          
          
          
          # the phylogenetic reconstruction with a slow one-state selection model
          transitionmatrix<-reconstruction_selection_four
          if(length(table(cladeFsimulatedtips))==3) transitionmatrix<-reconstruction_selection_three
          if(length(table(cladeFsimulatedtips))==2) transitionmatrix<-reconstruction_selection_two
          if(length(table(cladeFsimulatedtips))==1) transitionmatrix<-reconstruction_selection_one
          
          selection_to_a_reconstruction_twofold <- ace(cladeFsimulatedtips, cladeFtree, type="discrete", model=transitionmatrix)
          
          drift_results[counter,8]<-selection_to_a_reconstruction_twofold$loglik
          
          
          # the phylogenetic reconstruction with a pathway slow selection model
          transitionmatrix<-reconstruction_selection_four_pathway
          if(length(table(cladeFsimulatedtips))==3) transitionmatrix<-reconstruction_selection_three_pathway
          if(length(table(cladeFsimulatedtips))==2) transitionmatrix<-reconstruction_selection_two_pathway
          if(length(table(cladeFsimulatedtips))==1) transitionmatrix<-reconstruction_selection_one_pathway
          
          
          selection_to_a_reconstruction_pathway_slow <- fitDiscrete(dat=cladeFsimulatedtips, phy=cladeFtree, model="meristic")
          
          drift_results[counter,9]<-selection_to_a_reconstruction_pathway_slow$opt$lnL
          
          if(length(unique(cladeFsimulatedtips))>1) ifelse(names(table(cladeFsimulatedtips))[2]!=2, drift_results[counter,21]<-"TRUE",drift_results[counter,21]<-"FALSE")
          
          
          # the phylogenetic reconstruction with an alternative pathway selection model
          transitionmatrix<-reconstruction_selection_four_pathway_toone
          if(length(table(cladeFsimulatedtips))==3) transitionmatrix<-reconstruction_selection_three_pathway_toone
          if(length(table(cladeFsimulatedtips))==2) transitionmatrix<-reconstruction_selection_two_pathway_toone
          if(length(table(cladeFsimulatedtips))==1) transitionmatrix<-reconstruction_selection_one_pathway_toone
          
          selection_to_a_reconstruction_pathway_toone <- ace(cladeFsimulatedtips, cladeFtree, type="discrete", model=transitionmatrix)
          
          drift_results[counter,10]<-selection_to_a_reconstruction_pathway_toone$loglik
          
          # the phylogenetic reconstruction with an alternative pathway selection model
          transitionmatrix<-reconstruction_selection_four_pathway_fromone
          if(length(table(cladeFsimulatedtips))==3) transitionmatrix<-reconstruction_selection_three_pathway_fromone
          if(length(table(cladeFsimulatedtips))==2) transitionmatrix<-reconstruction_selection_two_pathway_fromone
          if(length(table(cladeFsimulatedtips))==1) transitionmatrix<-reconstruction_selection_one_pathway_fromone
          
          selection_to_a_reconstruction_pathway_fromone <- ace(cladeFsimulatedtips, cladeFtree, type="discrete", model=transitionmatrix)
          
          drift_results[counter,18]<-selection_to_a_reconstruction_pathway_fromone$loglik
          
          
          
          # the phylogenetic reconstruction with an unconstrained model
          unconstrained_reconstruction <- ace(cladeFsimulatedtips, cladeFtree, type="discrete", model="ARD")
          
          drift_results[counter,11]<-unconstrained_reconstruction$loglik
          
          resultphylosiglambda<-phylosig(cladeFtree,as.numeric(cladeFsimulatedtips),method="lambda") 
          drift_results[counter,12]<-resultphylosiglambda$lambda
          resultphylosigK<-phylosig(cladeFtree,as.numeric(cladeFsimulatedtips),method="K") 
          drift_results[counter,13]<-resultphylosigK[1]
          
        }
        
        counter<-counter+1
        
        
        
        
        
#--#--#--#--#--#--#        #-#-#-#-#-#-#-#-#-#-#-#
        
        # Subset the tree and the data to assume that 25% of data is randomly missing
        
        for (randomsample in 1:1) {
          
          ThreeQuarterClade <- currenttree$tip.label
          ThreeQuarterClade <- sample(ThreeQuarterClade)
          ThreeQuarterClade <- ThreeQuarterClade[1:129]
          
          ThreeQuartertree<-keep.tip(currenttree,ThreeQuarterClade)
          
          ThreeQuartersimulatedtips<-simulatedtips_neutral[names(simulatedtips_neutral) %in% ThreeQuarterClade]
          
          ThreeQuartersimulatedtips<-droplevels(ThreeQuartersimulatedtips)
          
          drift_results[counter,1]<-tree_used
          ifelse(drift_variant<4,drift_results[counter,2]<-2,drift_results[counter,2]<-4)
          ifelse(drift_variant<4,drift_results[counter,3]<-drift_variant,drift_results[counter,3]<-drift_variant-3)
          drift_results[counter,4]<-repetition
          
          drift_results[counter,5]<-"ThreeQuarterSample"
          drift_results[counter,6]<-length(unique(ThreeQuartersimulatedtips))
          
          if (length(table(ThreeQuartersimulatedtips))!=1) { 
            
            # the phylogenetic reconstruction with the drift model, the model under which the data were simulated
            equal_reconstruction <- ace(ThreeQuartersimulatedtips, ThreeQuartertree, type="discrete", model="ER")
            
            drift_results[counter,7]<-equal_reconstruction$loglik
            
            
            
            # the phylogenetic reconstruction with a slow one-state selection model
            transitionmatrix<-reconstruction_selection_four
            if(length(table(ThreeQuartersimulatedtips))==3) transitionmatrix<-reconstruction_selection_three
            if(length(table(ThreeQuartersimulatedtips))==2) transitionmatrix<-reconstruction_selection_two
            if(length(table(ThreeQuartersimulatedtips))==1) transitionmatrix<-reconstruction_selection_one
            
            selection_to_a_reconstruction_twofold <- ace(ThreeQuartersimulatedtips, ThreeQuartertree, type="discrete", model=transitionmatrix)
            
            drift_results[counter,8]<-selection_to_a_reconstruction_twofold$loglik
            
            
            # the phylogenetic reconstruction with a pathway slow selection model
            transitionmatrix<-reconstruction_selection_four_pathway
            if(length(table(ThreeQuartersimulatedtips))==3) transitionmatrix<-reconstruction_selection_three_pathway
            if(length(table(ThreeQuartersimulatedtips))==2) transitionmatrix<-reconstruction_selection_two_pathway
            if(length(table(ThreeQuartersimulatedtips))==1) transitionmatrix<-reconstruction_selection_one_pathway
            
            
            selection_to_a_reconstruction_pathway_slow <- fitDiscrete(dat=ThreeQuartersimulatedtips, phy=ThreeQuartertree, model="meristic")
            
            drift_results[counter,9]<-selection_to_a_reconstruction_pathway_slow$opt$lnL
            
            if(length(unique(ThreeQuartersimulatedtips))>1) ifelse(names(table(ThreeQuartersimulatedtips))[2]!=2, drift_results[counter,21]<-"TRUE",drift_results[counter,21]<-"FALSE")
            
            
            # the phylogenetic reconstruction with an alternative pathway selection model
            transitionmatrix<-reconstruction_selection_four_pathway_toone
            if(length(table(ThreeQuartersimulatedtips))==3) transitionmatrix<-reconstruction_selection_three_pathway_toone
            if(length(table(ThreeQuartersimulatedtips))==2) transitionmatrix<-reconstruction_selection_two_pathway_toone
            if(length(table(ThreeQuartersimulatedtips))==1) transitionmatrix<-reconstruction_selection_one_pathway_toone
            
            selection_to_a_reconstruction_pathway_toone <- ace(ThreeQuartersimulatedtips, ThreeQuartertree, type="discrete", model=transitionmatrix)
            
            drift_results[counter,10]<-selection_to_a_reconstruction_pathway_toone$loglik
            
            # the phylogenetic reconstruction with an alternative pathway selection model
            transitionmatrix<-reconstruction_selection_four_pathway_fromone
            if(length(table(ThreeQuartersimulatedtips))==3) transitionmatrix<-reconstruction_selection_three_pathway_fromone
            if(length(table(ThreeQuartersimulatedtips))==2) transitionmatrix<-reconstruction_selection_two_pathway_fromone
            if(length(table(ThreeQuartersimulatedtips))==1) transitionmatrix<-reconstruction_selection_one_pathway_fromone
            
            selection_to_a_reconstruction_pathway_fromone <- ace(ThreeQuartersimulatedtips, ThreeQuartertree, type="discrete", model=transitionmatrix)
            
            drift_results[counter,18]<-selection_to_a_reconstruction_pathway_fromone$loglik
            
            
            
            # the phylogenetic reconstruction with an unconstrained model
            unconstrained_reconstruction <- ace(ThreeQuartersimulatedtips, ThreeQuartertree, type="discrete", model="ARD")
            
            drift_results[counter,11]<-unconstrained_reconstruction$loglik
            
            resultphylosiglambda<-phylosig(ThreeQuartertree,as.numeric(ThreeQuartersimulatedtips),method="lambda") 
            drift_results[counter,12]<-resultphylosiglambda$lambda
            resultphylosigK<-phylosig(ThreeQuartertree,as.numeric(ThreeQuartersimulatedtips),method="K") 
            drift_results[counter,13]<-resultphylosigK[1]
            
          }
          
          counter<-counter+1
          
        } # end of the 10 random subsamples of 75% of societies
        
        
        
        
        
        #-#-#-#-#-#-#-#-#-#-#-#
        
        # Subset the tree and the data to assume that 50% of data is randomly missing
        
        for (randomsample in 1:1) {
          
          HalfClade <- currenttree$tip.label
          HalfClade <- sample(HalfClade)
          HalfClade <- HalfClade[1:86]
          
          Halftree<-keep.tip(currenttree,HalfClade)
          
          Halfsimulatedtips<-simulatedtips_neutral[names(simulatedtips_neutral) %in% HalfClade]
          
          Halfsimulatedtips<-droplevels(Halfsimulatedtips)
          
          drift_results[counter,1]<-tree_used
          ifelse(drift_variant<4,drift_results[counter,2]<-2,drift_results[counter,2]<-4)
          ifelse(drift_variant<4,drift_results[counter,3]<-drift_variant,drift_results[counter,3]<-drift_variant-3)
          drift_results[counter,4]<-repetition
          
          drift_results[counter,5]<-"HalfSample"
          drift_results[counter,6]<-length(unique(Halfsimulatedtips))
          
          if (length(table(Halfsimulatedtips))!=1) { 
            
            # the phylogenetic reconstruction with the drift model, the model under which the data were simulated
            equal_reconstruction <- ace(Halfsimulatedtips, Halftree, type="discrete", model="ER")
            
            drift_results[counter,7]<-equal_reconstruction$loglik
            
            
            
            # the phylogenetic reconstruction with a slow one-state selection model
            transitionmatrix<-reconstruction_selection_four
            if(length(table(Halfsimulatedtips))==3) transitionmatrix<-reconstruction_selection_three
            if(length(table(Halfsimulatedtips))==2) transitionmatrix<-reconstruction_selection_two
            if(length(table(Halfsimulatedtips))==1) transitionmatrix<-reconstruction_selection_one
            
            selection_to_a_reconstruction_twofold <- ace(Halfsimulatedtips, Halftree, type="discrete", model=transitionmatrix)
            
            drift_results[counter,8]<-selection_to_a_reconstruction_twofold$loglik
            
            
            # the phylogenetic reconstruction with a pathway slow selection model
            transitionmatrix<-reconstruction_selection_four_pathway
            if(length(table(Halfsimulatedtips))==3) transitionmatrix<-reconstruction_selection_three_pathway
            if(length(table(Halfsimulatedtips))==2) transitionmatrix<-reconstruction_selection_two_pathway
            if(length(table(Halfsimulatedtips))==1) transitionmatrix<-reconstruction_selection_one_pathway
            
            
            selection_to_a_reconstruction_pathway_slow <- fitDiscrete(dat=Halfsimulatedtips, phy=Halftree, model="meristic")
            
            drift_results[counter,9]<-selection_to_a_reconstruction_pathway_slow$opt$lnL
            
            
            if(length(unique(Halfsimulatedtips))>1) ifelse(names(table(Halfsimulatedtips))[2]!=2, drift_results[counter,21]<-"TRUE",drift_results[counter,21]<-"FALSE")
            
            
            # the phylogenetic reconstruction with an alternative pathway selection model
            transitionmatrix<-reconstruction_selection_four_pathway_toone
            if(length(table(Halfsimulatedtips))==3) transitionmatrix<-reconstruction_selection_three_pathway_toone
            if(length(table(Halfsimulatedtips))==2) transitionmatrix<-reconstruction_selection_two_pathway_toone
            if(length(table(Halfsimulatedtips))==1) transitionmatrix<-reconstruction_selection_one_pathway_toone
            
            selection_to_a_reconstruction_pathway_toone <- ace(Halfsimulatedtips, Halftree, type="discrete", model=transitionmatrix)
            
            drift_results[counter,10]<-selection_to_a_reconstruction_pathway_toone$loglik
            
            # the phylogenetic reconstruction with an alternative pathway selection model
            transitionmatrix<-reconstruction_selection_four_pathway_fromone
            if(length(table(Halfsimulatedtips))==3) transitionmatrix<-reconstruction_selection_three_pathway_fromone
            if(length(table(Halfsimulatedtips))==2) transitionmatrix<-reconstruction_selection_two_pathway_fromone
            if(length(table(Halfsimulatedtips))==1) transitionmatrix<-reconstruction_selection_one_pathway_fromone
            
            selection_to_a_reconstruction_pathway_fromone <- ace(Halfsimulatedtips, Halftree, type="discrete", model=transitionmatrix)
            
            drift_results[counter,18]<-selection_to_a_reconstruction_pathway_fromone$loglik
            
            
            
            # the phylogenetic reconstruction with an unconstrained model
            unconstrained_reconstruction <- ace(Halfsimulatedtips, Halftree, type="discrete", model="ARD")
            
            drift_results[counter,11]<-unconstrained_reconstruction$loglik
            
            resultphylosiglambda<-phylosig(Halftree,as.numeric(Halfsimulatedtips),method="lambda") 
            drift_results[counter,12]<-resultphylosiglambda$lambda
            resultphylosigK<-phylosig(Halftree,as.numeric(Halfsimulatedtips),method="K") 
            drift_results[counter,13]<-resultphylosigK[1]
            
          }
          
          counter<-counter+1
          
        } # end of the 10 random subsamples of 50% of societies
        
        
        
        
        #-#-#-#-#-#-#-#-#-#-#-#
        
        # State dependent loss of samples - 50% from population of least frequent variant
        
        for (randomsample in 1:1) {
          
          randomdeletion<-function(x) (x*rbinom(1,1,0.5))
          
          if(length(table(simulatedtips_neutral))==1) {
            onetips<-as.numeric(simulatedtips_neutral)
            names(onetips)<-names(simulatedtips_neutral)
            onetips<-sapply(as.numeric(onetips),randomdeletion)
            Rarelostsimulatedtips<-simulatedtips_neutral[onetips==0]
            Rarelost<-names(Rarelostsimulatedtips)
            Rarelosttree<-keep.tip(currenttree,Rarelost)
          }
          
          if(length(table(simulatedtips_neutral))==2) {
            twotips<-as.numeric(simulatedtips_neutral)
            names(twotips)<-names(simulatedtips_neutral)
            twotips[twotips %in% c(1)]<-0
            twotips[twotips %in% c(2)]<-1
            twotips<-sapply(as.numeric(twotips),randomdeletion)
            Rarelostsimulatedtips<-simulatedtips_neutral[twotips==0]
            Rarelost<-names(Rarelostsimulatedtips)
            Rarelosttree<-keep.tip(currenttree,Rarelost)
          }
          
          if(length(table(simulatedtips_neutral))==3) {
            threetips<-as.numeric(simulatedtips_neutral)
            names(threetips)<-names(simulatedtips_neutral)
            threetips[threetips %in% c(1,2)]<-0
            threetips[threetips %in% c(3)]<-1
            threetips<-sapply(as.numeric(threetips),randomdeletion)
            Rarelostsimulatedtips<-simulatedtips_neutral[threetips==0]
            Rarelost<-names(Rarelostsimulatedtips)
            Rarelosttree<-keep.tip(currenttree,Rarelost)
          }    
          
          if(length(table(simulatedtips_neutral))==4) {
            fourtips<-as.numeric(simulatedtips_neutral)
            names(fourtips)<-names(simulatedtips_neutral)
            fourtips[fourtips %in% c(1,2,3)]<-0
            fourtips[fourtips %in% c(4)]<-1
            fourtips<-sapply(as.numeric(fourtips),randomdeletion)
            Rarelostsimulatedtips<-simulatedtips_neutral[fourtips==0]
            Rarelost<-names(Rarelostsimulatedtips)
            Rarelosttree<-keep.tip(currenttree,Rarelost)
          }    
          
          Rarelostsimulatedtips<-droplevels(Rarelostsimulatedtips)
          
          drift_results[counter,1]<-tree_used
          ifelse(drift_variant<4,drift_results[counter,2]<-2,drift_results[counter,2]<-4)
          ifelse(drift_variant<4,drift_results[counter,3]<-drift_variant,drift_results[counter,3]<-drift_variant-3)
          drift_results[counter,4]<-repetition
          
          drift_results[counter,5]<-"Rarelost"
          drift_results[counter,6]<-length(unique(Rarelostsimulatedtips))
          
          if (length(table(Rarelostsimulatedtips))!=1) { 
            
            # the phylogenetic reconstruction with the drift model, the model under which the data were simulated
            equal_reconstruction <- ace(Rarelostsimulatedtips, Rarelosttree, type="discrete", model="ER")
            
            drift_results[counter,7]<-equal_reconstruction$loglik
            
            
            
            # the phylogenetic reconstruction with a slow one-state selection model
            transitionmatrix<-reconstruction_selection_four
            if(length(table(Rarelostsimulatedtips))==3) transitionmatrix<-reconstruction_selection_three
            if(length(table(Rarelostsimulatedtips))==2) transitionmatrix<-reconstruction_selection_two
            if(length(table(Rarelostsimulatedtips))==1) transitionmatrix<-reconstruction_selection_one
            
            selection_to_a_reconstruction_twofold <- ace(Rarelostsimulatedtips, Rarelosttree, type="discrete", model=transitionmatrix)
            
            drift_results[counter,8]<-selection_to_a_reconstruction_twofold$loglik
            
            
            # the phylogenetic reconstruction with a pathway slow selection model
            transitionmatrix<-reconstruction_selection_four_pathway
            if(length(table(Rarelostsimulatedtips))==3) transitionmatrix<-reconstruction_selection_three_pathway
            if(length(table(Rarelostsimulatedtips))==2) transitionmatrix<-reconstruction_selection_two_pathway
            if(length(table(Rarelostsimulatedtips))==1) transitionmatrix<-reconstruction_selection_one_pathway
            
            selection_to_a_reconstruction_pathway_slow <- fitDiscrete(dat=Rarelostsimulatedtips, phy=Rarelosttree, model="meristic")
            
            drift_results[counter,9]<-selection_to_a_reconstruction_pathway_slow$opt$lnL
            
            if(length(unique(Rarelostsimulatedtips))>1) ifelse(names(table(Rarelostsimulatedtips))[2]!=2, drift_results[counter,21]<-"TRUE",drift_results[counter,21]<-"FALSE")
            
            
            # the phylogenetic reconstruction with an alternative pathway selection model
            transitionmatrix<-reconstruction_selection_four_pathway_toone
            if(length(table(Rarelostsimulatedtips))==3) transitionmatrix<-reconstruction_selection_three_pathway_toone
            if(length(table(Rarelostsimulatedtips))==2) transitionmatrix<-reconstruction_selection_two_pathway_toone
            if(length(table(Rarelostsimulatedtips))==1) transitionmatrix<-reconstruction_selection_one_pathway_toone
            
            selection_to_a_reconstruction_pathway_toone <- ace(Rarelostsimulatedtips, Rarelosttree, type="discrete", model=transitionmatrix)
            
            drift_results[counter,10]<-selection_to_a_reconstruction_pathway_toone$loglik
            
            # the phylogenetic reconstruction with an alternative pathway selection model
            transitionmatrix<-reconstruction_selection_four_pathway_fromone
            if(length(table(Rarelostsimulatedtips))==3) transitionmatrix<-reconstruction_selection_three_pathway_fromone
            if(length(table(Rarelostsimulatedtips))==2) transitionmatrix<-reconstruction_selection_two_pathway_fromone
            if(length(table(Rarelostsimulatedtips))==1) transitionmatrix<-reconstruction_selection_one_pathway_fromone
            
            selection_to_a_reconstruction_pathway_fromone <- ace(Rarelostsimulatedtips, Rarelosttree, type="discrete", model=transitionmatrix)
            
            drift_results[counter,18]<-selection_to_a_reconstruction_pathway_fromone$loglik
            
            
            
            # the phylogenetic reconstruction with an unconstrained model
            unconstrained_reconstruction <- ace(Rarelostsimulatedtips, Rarelosttree, type="discrete", model="ARD")
            
            drift_results[counter,11]<-unconstrained_reconstruction$loglik
            
            
            resultphylosiglambda<-phylosig(Rarelosttree,as.numeric(Rarelostsimulatedtips),method="lambda") 
            drift_results[counter,12]<-resultphylosiglambda$lambda
            resultphylosigK<-phylosig(Rarelosttree,as.numeric(Rarelostsimulatedtips),method="K") 
            drift_results[counter,13]<-resultphylosigK[1]
            
          }
          
          counter<-counter+1
          
        } # end of the 10 random loss of 50% of societies with the rare variant
        
        
        
        
        #-#-#-#-#-#-#-#-#-#-#-#
        
        # State dependent loss of samples - 50% from population of most frequent variant
        
        for (randomsample in 1:1) {
          
          randomdeletion<-function(x) (x*rbinom(1,1,0.5))
          
          if(length(table(simulatedtips_neutral))==1) {
            onetips<-as.numeric(simulatedtips_neutral)
            names(onetips)<-names(simulatedtips_neutral)
            onetips<-sapply(as.numeric(onetips),randomdeletion)
            Frequentlostsimulatedtips<-simulatedtips_neutral[onetips==0]
            Frequentlost<-names(Frequentlostsimulatedtips)
            Frequentlosttree<-keep.tip(currenttree,Frequentlost)
          }
          
          if(length(table(simulatedtips_neutral))==2) {
            twotips<-as.numeric(simulatedtips_neutral)
            names(twotips)<-names(simulatedtips_neutral)
            twotips[twotips %in% c(2)]<-0
            twotips[twotips %in% c(1)]<-1
            twotips<-sapply(as.numeric(twotips),randomdeletion)
            Frequentlostsimulatedtips<-simulatedtips_neutral[twotips==0]
            Frequentlost<-names(Frequentlostsimulatedtips)
            Frequentlosttree<-keep.tip(currenttree,Frequentlost)
          }
          
          if(length(table(simulatedtips_neutral))==3) {
            threetips<-as.numeric(simulatedtips_neutral)
            names(threetips)<-names(simulatedtips_neutral)
            threetips[threetips %in% c(3,2)]<-0
            threetips[threetips %in% c(1)]<-1
            threetips<-sapply(as.numeric(threetips),randomdeletion)
            Frequentlostsimulatedtips<-simulatedtips_neutral[threetips==0]
            Frequentlost<-names(Frequentlostsimulatedtips)
            Frequentlosttree<-keep.tip(currenttree,Frequentlost)
          }    
          
          if(length(table(simulatedtips_neutral))==4) {
            fourtips<-as.numeric(simulatedtips_neutral)
            names(fourtips)<-names(simulatedtips_neutral)
            fourtips[fourtips %in% c(4,2,3)]<-0
            fourtips[fourtips %in% c(1)]<-1
            fourtips<-sapply(as.numeric(fourtips),randomdeletion)
            Frequentlostsimulatedtips<-simulatedtips_neutral[fourtips==0]
            Frequentlost<-names(Frequentlostsimulatedtips)
            Frequentlosttree<-keep.tip(currenttree,Frequentlost)
          }    
          
          Frequentlostsimulatedtips<-droplevels(Frequentlostsimulatedtips)
          
          drift_results[counter,1]<-tree_used
          ifelse(drift_variant<4,drift_results[counter,2]<-2,drift_results[counter,2]<-4)
          ifelse(drift_variant<4,drift_results[counter,3]<-drift_variant,drift_results[counter,3]<-drift_variant-3)
          drift_results[counter,4]<-repetition
          
          drift_results[counter,5]<-"Frequentlost"
          drift_results[counter,6]<-length(unique(Frequentlostsimulatedtips))
          
          if (length(table(Frequentlostsimulatedtips))!=1) { 
            
            # the phylogenetic reconstruction with the drift model, the model under which the data were simulated
            equal_reconstruction <- ace(Frequentlostsimulatedtips, Frequentlosttree, type="discrete", model="ER")
            
            drift_results[counter,7]<-equal_reconstruction$loglik
            
            
            
            # the phylogenetic reconstruction with a slow one-state selection model
            transitionmatrix<-reconstruction_selection_four
            if(length(table(Frequentlostsimulatedtips))==3) transitionmatrix<-reconstruction_selection_three
            if(length(table(Frequentlostsimulatedtips))==2) transitionmatrix<-reconstruction_selection_two
            if(length(table(Frequentlostsimulatedtips))==1) transitionmatrix<-reconstruction_selection_one
            
            selection_to_a_reconstruction_twofold <- ace(Frequentlostsimulatedtips, Frequentlosttree, type="discrete", model=transitionmatrix)
            
            drift_results[counter,8]<-selection_to_a_reconstruction_twofold$loglik
            
            
            # the phylogenetic reconstruction with a pathway slow selection model
            transitionmatrix<-reconstruction_selection_four_pathway
            if(length(table(Frequentlostsimulatedtips))==3) transitionmatrix<-reconstruction_selection_three_pathway
            if(length(table(Frequentlostsimulatedtips))==2) transitionmatrix<-reconstruction_selection_two_pathway
            if(length(table(Frequentlostsimulatedtips))==1) transitionmatrix<-reconstruction_selection_one_pathway
            
            selection_to_a_reconstruction_pathway_slow <- fitDiscrete(dat=Frequentlostsimulatedtips, phy=Frequentlosttree, model="meristic")
            
            drift_results[counter,9]<-selection_to_a_reconstruction_pathway_slow$opt$lnL
            
            if(length(unique(Frequentlostsimulatedtips))>1) ifelse(names(table(Frequentlostsimulatedtips))[2]!=2, drift_results[counter,21]<-"TRUE",drift_results[counter,21]<-"FALSE")
            
            
            # the phylogenetic reconstruction with an alternative pathway selection model
            transitionmatrix<-reconstruction_selection_four_pathway_toone
            if(length(table(Frequentlostsimulatedtips))==3) transitionmatrix<-reconstruction_selection_three_pathway_toone
            if(length(table(Frequentlostsimulatedtips))==2) transitionmatrix<-reconstruction_selection_two_pathway_toone
            if(length(table(Frequentlostsimulatedtips))==1) transitionmatrix<-reconstruction_selection_one_pathway_toone
            
            selection_to_a_reconstruction_pathway_toone <- ace(Frequentlostsimulatedtips, Frequentlosttree, type="discrete", model=transitionmatrix)
            
            drift_results[counter,10]<-selection_to_a_reconstruction_pathway_toone$loglik
            
            # the phylogenetic reconstruction with an alternative pathway selection model
            transitionmatrix<-reconstruction_selection_four_pathway_fromone
            if(length(table(Frequentlostsimulatedtips))==3) transitionmatrix<-reconstruction_selection_three_pathway_fromone
            if(length(table(Frequentlostsimulatedtips))==2) transitionmatrix<-reconstruction_selection_two_pathway_fromone
            if(length(table(Frequentlostsimulatedtips))==1) transitionmatrix<-reconstruction_selection_one_pathway_fromone
            
            selection_to_a_reconstruction_pathway_fromone <- ace(Frequentlostsimulatedtips, Frequentlosttree, type="discrete", model=transitionmatrix)
            
            drift_results[counter,18]<-selection_to_a_reconstruction_pathway_fromone$loglik
            
            
            
            # the phylogenetic reconstruction with an unconstrained model
            unconstrained_reconstruction <- ace(Frequentlostsimulatedtips, Frequentlosttree, type="discrete", model="ARD")
            
            drift_results[counter,11]<-unconstrained_reconstruction$loglik
            
            
            resultphylosiglambda<-phylosig(Frequentlosttree,as.numeric(Frequentlostsimulatedtips),method="lambda") 
            drift_results[counter,12]<-resultphylosiglambda$lambda
            resultphylosigK<-phylosig(Frequentlosttree,as.numeric(Frequentlostsimulatedtips),method="K") 
            drift_results[counter,13]<-resultphylosigK[1]
            
          }
          
          counter<-counter+1
          
        } # end of the 10 random loss of 50% of societies with the frequent variant
        
        
        
        
        
        # Loss according to ecoregion (50% of those living in forests)
        
        
        # Subset the tree and the data to assume that 50% of data is randomly missing
        
        
        for (randomsample in 1:1) {
          
          randomdeletion<-function(x) (x*rbinom(1,1,0.5))
          Foresttips<-sapply(as.numeric(Forestmember),randomdeletion)
          Forestlosstips<-simulatedtips_neutral[Foresttips==0]
          Forestlosttree<-keep.tip(currenttree,names(Forestlosstips))
          
          Forestlosstips<-droplevels(Forestlosstips)
          
          drift_results[counter,1]<-tree_used
          ifelse(drift_variant<4,drift_results[counter,2]<-2,drift_results[counter,2]<-4)
          ifelse(drift_variant<4,drift_results[counter,3]<-drift_variant,drift_results[counter,3]<-drift_variant-3)
          drift_results[counter,4]<-repetition
          
          drift_results[counter,5]<-"Forestloss"
          drift_results[counter,6]<-length(unique(Forestlosstips))
          
          if (length(table(Forestlosstips))!=1) { 
            
            # the phylogenetic reconstruction with the drift model, the model under which the data were simulated
            equal_reconstruction <- ace(Forestlosstips, Forestlosttree, type="discrete", model="ER")
            
            drift_results[counter,7]<-equal_reconstruction$loglik
            
            
            
            # the phylogenetic reconstruction with a slow one-state selection model
            transitionmatrix<-reconstruction_selection_four
            if(length(table(Forestlosstips))==3) transitionmatrix<-reconstruction_selection_three
            if(length(table(Forestlosstips))==2) transitionmatrix<-reconstruction_selection_two
            if(length(table(Forestlosstips))==1) transitionmatrix<-reconstruction_selection_one
            
            selection_to_a_reconstruction_twofold <- ace(Forestlosstips, Forestlosttree, type="discrete", model=transitionmatrix)
            
            drift_results[counter,8]<-selection_to_a_reconstruction_twofold$loglik
            
            
            # the phylogenetic reconstruction with a pathway slow selection model
            transitionmatrix<-reconstruction_selection_four_pathway
            if(length(table(Forestlosstips))==3) transitionmatrix<-reconstruction_selection_three_pathway
            if(length(table(Forestlosstips))==2) transitionmatrix<-reconstruction_selection_two_pathway
            if(length(table(Forestlosstips))==1) transitionmatrix<-reconstruction_selection_one_pathway
            
            selection_to_a_reconstruction_pathway_slow <- fitDiscrete(dat=Forestlosstips, phy=Forestlosttree, model="meristic")
            
            drift_results[counter,9]<-selection_to_a_reconstruction_pathway_slow$opt$lnL
            
            if(length(unique(Forestlosstips))>1) ifelse(names(table(Forestlosstips))[2]!=2, drift_results[counter,21]<-"TRUE",drift_results[counter,21]<-"FALSE") 
            
            
            # the phylogenetic reconstruction with an alternative pathway selection model
            transitionmatrix<-reconstruction_selection_four_pathway_toone
            if(length(table(Forestlosstips))==3) transitionmatrix<-reconstruction_selection_three_pathway_toone
            if(length(table(Forestlosstips))==2) transitionmatrix<-reconstruction_selection_two_pathway_toone
            if(length(table(Forestlosstips))==1) transitionmatrix<-reconstruction_selection_one_pathway_toone
            
            selection_to_a_reconstruction_pathway_toone <- ace(Forestlosstips, Forestlosttree, type="discrete", model=transitionmatrix)
            
            drift_results[counter,10]<-selection_to_a_reconstruction_pathway_toone$loglik
            
            # the phylogenetic reconstruction with an alternative pathway selection model
            transitionmatrix<-reconstruction_selection_four_pathway_fromone
            if(length(table(Forestlosstips))==3) transitionmatrix<-reconstruction_selection_three_pathway_fromone
            if(length(table(Forestlosstips))==2) transitionmatrix<-reconstruction_selection_two_pathway_fromone
            if(length(table(Forestlosstips))==1) transitionmatrix<-reconstruction_selection_one_pathway_fromone
            
            selection_to_a_reconstruction_pathway_fromone <- ace(Forestlosstips, Forestlosttree, type="discrete", model=transitionmatrix)
            
            drift_results[counter,18]<-selection_to_a_reconstruction_pathway_fromone$loglik
            
            
            
            # the phylogenetic reconstruction with an unconstrained model
            unconstrained_reconstruction <- ace(Forestlosstips, Forestlosttree, type="discrete", model="ARD")
            
            drift_results[counter,11]<-unconstrained_reconstruction$loglik
            
            resultphylosiglambda<-phylosig(Forestlosttree,as.numeric(Forestlosstips),method="lambda") 
            drift_results[counter,12]<-resultphylosiglambda$lambda
            resultphylosigK<-phylosig(Forestlosttree,as.numeric(Forestlosstips),method="K") 
            drift_results[counter,13]<-resultphylosigK[1]
            
          }
          
          counter<-counter+1
          
        } # end of the 10 random subsamples of 50% of societies lost from forests
        
        

        
        
                
#--#--#--#--#--#--#        #-#-#-#-#-#-#-#-#-#-#-#-#
        
        
        # Horizontal transfer, 10% of societies change, randomly adopt another state
        
        for (randomsample in 1:1) {
          
          
          societiestochange<-sample(simulatedtips_neutral,17)
          societiestochange[societiestochange==1]<-NA
          societiestochange[societiestochange==2]<-1
          societiestochange[is.na(societiestochange)==T]<-2
          horizontal10tips<-simulatedtips_neutral
          horizontal10tips[names(horizontal10tips) %in% names(societiestochange)]<-societiestochange
          
          horizontal10tips<-droplevels(horizontal10tips)
          
          
          #Start filling in the data in the respective column
          drift_results[counter,1]<-tree_used
          ifelse(drift_variant<4,drift_results[counter,2]<-2,drift_results[counter,2]<-4)
          ifelse(drift_variant<4,drift_results[counter,3]<-drift_variant,drift_results[counter,3]<-drift_variant-3)
          drift_results[counter,4]<-repetition
          
          
          #Start with the full sample of 172 societies
          
          drift_results[counter,5]<-"horizontal10"
          drift_results[counter,6]<-length(unique(horizontal10tips))
          
          if (length(table(horizontal10tips))!=1) { 
            
            # the phylogenetic reconstruction with the drift model, the model under which the data were simulated
            equal_reconstruction <- ace(horizontal10tips, currenttree, type="discrete", model="ER")
            
            drift_results[counter,7]<-equal_reconstruction$loglik
            
            
            
            # the phylogenetic reconstruction with a one-state selection model
            transitionmatrix<-reconstruction_selection_four
            if(length(table(horizontal10tips))==3) transitionmatrix<-reconstruction_selection_three
            if(length(table(horizontal10tips))==2) transitionmatrix<-reconstruction_selection_two
            if(length(table(horizontal10tips))==1) transitionmatrix<-reconstruction_selection_one
            
            selection_to_a_reconstruction_twofold <- ace(horizontal10tips, currenttree, type="discrete", model=transitionmatrix)
            
            drift_results[counter,8]<-selection_to_a_reconstruction_twofold$loglik
            
            
            
            # the phylogenetic reconstruction with a pathway selection model
            transitionmatrix<-reconstruction_selection_four_pathway
            if(length(table(horizontal10tips))==3) transitionmatrix<-reconstruction_selection_three_pathway
            if(length(table(horizontal10tips))==2) transitionmatrix<-reconstruction_selection_two_pathway
            if(length(table(horizontal10tips))==1) transitionmatrix<-reconstruction_selection_one_pathway
            
            
            selection_to_a_reconstruction_pathway_slow <- fitDiscrete(dat=horizontal10tips, phy=currenttree, model="meristic")
            
            drift_results[counter,9]<-selection_to_a_reconstruction_pathway_slow$opt$lnL
            
            if(length(unique(horizontal10tips))>1) ifelse(names(table(horizontal10tips))[2]!=2, drift_results[counter,21]<-"TRUE",drift_results[counter,21]<-"FALSE") 
            
            
            # the phylogenetic reconstruction with an alternative pathway selection model
            transitionmatrix<-reconstruction_selection_four_pathway_toone
            if(length(table(horizontal10tips))==3) transitionmatrix<-reconstruction_selection_three_pathway_toone
            if(length(table(horizontal10tips))==2) transitionmatrix<-reconstruction_selection_two_pathway_toone
            if(length(table(horizontal10tips))==1) transitionmatrix<-reconstruction_selection_one_pathway_toone
            
            selection_to_a_reconstruction_pathway_toone <- ace(horizontal10tips, currenttree, type="discrete", model=transitionmatrix)
            
            drift_results[counter,10]<-selection_to_a_reconstruction_pathway_toone$loglik
            
            # the phylogenetic reconstruction with an alternative pathway selection model
            transitionmatrix<-reconstruction_selection_four_pathway_fromone
            if(length(table(horizontal10tips))==3) transitionmatrix<-reconstruction_selection_three_pathway_fromone
            if(length(table(horizontal10tips))==2) transitionmatrix<-reconstruction_selection_two_pathway_fromone
            if(length(table(horizontal10tips))==1) transitionmatrix<-reconstruction_selection_one_pathway_fromone
            
            selection_to_a_reconstruction_pathway_fromone <- ace(horizontal10tips, currenttree, type="discrete", model=transitionmatrix)
            
            drift_results[counter,18]<-selection_to_a_reconstruction_pathway_fromone$loglik
            
            
            
            # the phylogenetic reconstruction with an unconstrained model
            unconstrained_reconstruction <- ace(horizontal10tips, currenttree, type="discrete", model="ARD")
            
            drift_results[counter,11]<-unconstrained_reconstruction$loglik
            
            resultphylosiglambda<-phylosig(currenttree,as.numeric(horizontal10tips),method="lambda") 
            drift_results[counter,12]<-resultphylosiglambda$lambda
            resultphylosigK<-phylosig(currenttree,as.numeric(horizontal10tips),method="K") 
            drift_results[counter,13]<-resultphylosigK[1]
            
            
            #For the simulations where there are only two variants, we check whether by chance the variants
            #end up isolated in Clade A or associated with societies living in forest ecoregions
            #For this we use Bayestraits Discrete, assessing the likelihood of independent versus dependent evolution
            if(length(table(horizontal10tips))==2) {
              
              bayestraitsdiscretedata<-matrix(NA,nrow=length(horizontal10tips),ncol=3)
              colnames(bayestraitsdiscretedata)<-c("Species","Simulated","CladeA")
              rownames(bayestraitsdiscretedata)<-c(1:length(horizontal10tips))
              bayestraitsdiscretedata[,1]<-names(horizontal10tips)
              bayestraitsdiscretedata[,2]<-as.numeric(as.character(horizontal10tips))
              bayestraitsdiscretedata[,2]<-as.numeric(bayestraitsdiscretedata[,2])-as.numeric(min(bayestraitsdiscretedata[,2]))
              bayestraitsdiscretedata[bayestraitsdiscretedata[,2]!=0,2]<-1
              bayestraitsdiscretedata[,3]<- as.numeric(as.character(cladeAmember))
              bayestraitsdiscretedata<-as.data.frame(bayestraitsdiscretedata)
              
              command_vec1 <- c("2", "1") #option 1 = 2 discrete independent; option 2 = 1 maximum likelihood
              results_1 <- bayestraits(bayestraitsdiscretedata, currenttree, command_vec1)
              log_1 <- results_1$Log
              drift_results[counter,14]<-log_1$results$Lh
              
              command_vec2 <- c("3", "1") #option 1 = 3 discrete dependent; option 2 = 1 maximum likelihood
              results_2 <- bayestraits(bayestraitsdiscretedata, currenttree, command_vec2)
              log_2 <- results_2$Log
              drift_results[counter,15]<-log_2$results$Lh
              
              
              bayestraitsdiscretedata<-matrix(NA,nrow=length(horizontal10tips),ncol=3)
              colnames(bayestraitsdiscretedata)<-c("Species","Simulated","Forestmember")
              rownames(bayestraitsdiscretedata)<-c(1:length(horizontal10tips))
              bayestraitsdiscretedata[,1]<-names(horizontal10tips)
              bayestraitsdiscretedata[,2]<-as.numeric(as.character(horizontal10tips))
              bayestraitsdiscretedata[,2]<-as.numeric(bayestraitsdiscretedata[,2])-as.numeric(min(bayestraitsdiscretedata[,2]))
              bayestraitsdiscretedata[bayestraitsdiscretedata[,2]!=0,2]<-1
              bayestraitsdiscretedata[,3]<- as.numeric(as.character(Forestmember))
              bayestraitsdiscretedata<-as.data.frame(bayestraitsdiscretedata)
              
              command_vec3 <- c("2", "1") #option 1 = 2 discrete independent; option 2 = 1 maximum likelihood
              results_3 <- bayestraits(bayestraitsdiscretedata, currenttree, command_vec3)
              log_3 <- results_3$Log
              drift_results[counter,16]<-log_3$results$Lh
              
              command_vec4 <- c("3", "1") #option 1 = 3 discrete dependent; option 2 = 1 maximum likelihood
              results_4 <- bayestraits(bayestraitsdiscretedata, currenttree, command_vec4)
              log_4 <- results_4$Log
              drift_results[counter,17]<-log_4$results$Lh
              
            } #end of bayestraits option
            
            
            #For the simulations where there are more than two variants, we run the Bayestraits multistate option
            
            if(length(table(horizontal10tips))>2) {
              
              bayestraitsdiscretedata<-matrix(NA,nrow=length(horizontal10tips),ncol=2)
              colnames(bayestraitsdiscretedata)<-c("Species","Simulated")
              rownames(bayestraitsdiscretedata)<-c(1:length(horizontal10tips))
              bayestraitsdiscretedata[,1]<-names(horizontal10tips)
              bayestraitsdiscretedata[,2]<-as.numeric(as.character(horizontal10tips))
              bayestraitsdiscretedata[,2]<-as.numeric(bayestraitsdiscretedata[,2])-1
              bayestraitsdiscretedata<-as.data.frame(bayestraitsdiscretedata)
              
              
              command_vec1 <- c("1", "1") #option 1 = 1 multistate; option 2 = 1 maximum likelihood
              results_1 <- bayestraits(bayestraitsdiscretedata, currenttree, command_vec1)
              log_1 <- results_1$Log
              drift_results[counter,19]<-log_1$results$Lh
              
              command_vec2 <- c("1", "1", "Restrict q01 q02 q03 q10 q12 q13 q20 q21 q23 q30 q31 q32") #option 1 = 1 multistate; option 2 = 1 maximum likelihood; option restrict all dependent
              results_2 <- bayestraits(bayestraitsdiscretedata, currenttree, command_vec2)
              log_2 <- results_2$Log
              drift_results[counter,20]<-log_2$results$Lh
              
              
            } #end of bayestraits multistate option
            
            
          }
          
          counter<-counter+1
          
        } # end of the 10 random subsamples with 10% horizontal transmission
        
        
        
        
        
        
        # Horizontal transfer, 10% of societies change, adopt state most frequent within their own clade
        
        for (randomsample in 1:1) {
          
          if (length(table(cladeAsimulatedtips))!=1) {
            
            cladeAtochange<-sample(cladeAsimulatedtips[cladeAsimulatedtips==names(table(cladeAsimulatedtips)[2])],min(17,table(cladeAsimulatedtips)[2]))
            levels(cladeAtochange)<-c("1","2","3","4")
            cladeAtochange[]<-1
            horizontalClade10tips<-simulatedtips_neutral
            horizontalClade10tips[names(horizontalClade10tips) %in% names(cladeAtochange)]<-cladeAtochange
            
            horizontalClade10tips<-droplevels(horizontalClade10tips)
            
            #Start filling in the data in the respective column
            drift_results[counter,1]<-tree_used
            ifelse(drift_variant<4,drift_results[counter,2]<-2,drift_results[counter,2]<-4)
            ifelse(drift_variant<4,drift_results[counter,3]<-drift_variant,drift_results[counter,3]<-drift_variant-3)
            drift_results[counter,4]<-repetition
            
            
            #Start with the full sample of 172 societies
            
            drift_results[counter,5]<-"horizontalClade10"
            drift_results[counter,6]<-length(unique(horizontalClade10tips))
            
            if (length(table(horizontalClade10tips))!=1) { 
              
              
              # the phylogenetic reconstruction with the drift model, the model under which the data were simulated
              equal_reconstruction <- ace(horizontalClade10tips, currenttree, type="discrete", model="ER")
              
              drift_results[counter,7]<-equal_reconstruction$loglik
              
              
              
              # the phylogenetic reconstruction with a one-state selection model
              transitionmatrix<-reconstruction_selection_four
              if(length(table(horizontalClade10tips))==3) transitionmatrix<-reconstruction_selection_three
              if(length(table(horizontalClade10tips))==2) transitionmatrix<-reconstruction_selection_two
              if(length(table(horizontalClade10tips))==1) transitionmatrix<-reconstruction_selection_one
              
              selection_to_a_reconstruction_twofold <- ace(horizontalClade10tips, currenttree, type="discrete", model=transitionmatrix)
              
              drift_results[counter,8]<-selection_to_a_reconstruction_twofold$loglik
              
              
              
              # the phylogenetic reconstruction with a pathway selection model
              transitionmatrix<-reconstruction_selection_four_pathway
              if(length(table(horizontalClade10tips))==3) transitionmatrix<-reconstruction_selection_three_pathway
              if(length(table(horizontalClade10tips))==2) transitionmatrix<-reconstruction_selection_two_pathway
              if(length(table(horizontalClade10tips))==1) transitionmatrix<-reconstruction_selection_one_pathway
              
              selection_to_a_reconstruction_pathway_slow <- fitDiscrete(dat=horizontalClade10tips, phy=currenttree, model="meristic")
              
              drift_results[counter,9]<-selection_to_a_reconstruction_pathway_slow$opt$lnL
              
              if(length(unique(horizontalClade10tips))>1) ifelse(names(table(horizontalClade10tips))[2]!=2, drift_results[counter,21]<-"TRUE",drift_results[counter,21]<-"FALSE") 
              
              
              
              
              # the phylogenetic reconstruction with an alternative pathway selection model
              transitionmatrix<-reconstruction_selection_four_pathway_toone
              if(length(table(horizontalClade10tips))==3) transitionmatrix<-reconstruction_selection_three_pathway_toone
              if(length(table(horizontalClade10tips))==2) transitionmatrix<-reconstruction_selection_two_pathway_toone
              if(length(table(horizontalClade10tips))==1) transitionmatrix<-reconstruction_selection_one_pathway_toone
              
              selection_to_a_reconstruction_pathway_toone <- ace(horizontalClade10tips, currenttree, type="discrete", model=transitionmatrix)
              
              drift_results[counter,10]<-selection_to_a_reconstruction_pathway_toone$loglik
              
              # the phylogenetic reconstruction with an alternative pathway selection model
              transitionmatrix<-reconstruction_selection_four_pathway_fromone
              if(length(table(horizontalClade10tips))==3) transitionmatrix<-reconstruction_selection_three_pathway_fromone
              if(length(table(horizontalClade10tips))==2) transitionmatrix<-reconstruction_selection_two_pathway_fromone
              if(length(table(horizontalClade10tips))==1) transitionmatrix<-reconstruction_selection_one_pathway_fromone
              
              selection_to_a_reconstruction_pathway_fromone <- ace(horizontalClade10tips, currenttree, type="discrete", model=transitionmatrix)
              
              drift_results[counter,18]<-selection_to_a_reconstruction_pathway_fromone$loglik
              
              
              
              # the phylogenetic reconstruction with an unconstrained model
              unconstrained_reconstruction <- ace(horizontalClade10tips, currenttree, type="discrete", model="ARD")
              
              drift_results[counter,11]<-unconstrained_reconstruction$loglik
              
              resultphylosiglambda<-phylosig(currenttree,as.numeric(horizontalClade10tips),method="lambda") 
              drift_results[counter,12]<-resultphylosiglambda$lambda
              resultphylosigK<-phylosig(currenttree,as.numeric(horizontalClade10tips),method="K") 
              drift_results[counter,13]<-resultphylosigK[1]
              
              
              #For the simulations where there are only two variants, we check whether by chance the variants
              #end up isolated in Clade A or associated with societies living in forest ecoregions
              #For this we use Bayestraits Discrete, assessing the likelihood of independent versus dependent evolution
              if(length(table(horizontalClade10tips))==2) {
                
                bayestraitsdiscretedata<-matrix(NA,nrow=length(horizontalClade10tips),ncol=3)
                colnames(bayestraitsdiscretedata)<-c("Species","Simulated","CladeA")
                rownames(bayestraitsdiscretedata)<-c(1:length(horizontalClade10tips))
                bayestraitsdiscretedata[,1]<-names(horizontalClade10tips)
                bayestraitsdiscretedata[,2]<-as.numeric(as.character(horizontalClade10tips))
                bayestraitsdiscretedata[,2]<-as.numeric(bayestraitsdiscretedata[,2])-as.numeric(min(bayestraitsdiscretedata[,2]))
                bayestraitsdiscretedata[bayestraitsdiscretedata[,2]!=0,2]<-1
                bayestraitsdiscretedata[,3]<- as.numeric(as.character(cladeAmember))
                bayestraitsdiscretedata<-as.data.frame(bayestraitsdiscretedata)
                
                command_vec1 <- c("2", "1") #option 1 = 2 discrete independent; option 2 = 1 maximum likelihood
                results_1 <- bayestraits(bayestraitsdiscretedata, currenttree, command_vec1)
                log_1 <- results_1$Log
                drift_results[counter,14]<-log_1$results$Lh
                
                command_vec2 <- c("3", "1") #option 1 = 3 discrete dependent; option 2 = 1 maximum likelihood
                results_2 <- bayestraits(bayestraitsdiscretedata, currenttree, command_vec2)
                log_2 <- results_2$Log
                drift_results[counter,15]<-log_2$results$Lh
                
                
                bayestraitsdiscretedata<-matrix(NA,nrow=length(horizontalClade10tips),ncol=3)
                colnames(bayestraitsdiscretedata)<-c("Species","Simulated","Forestmember")
                rownames(bayestraitsdiscretedata)<-c(1:length(horizontalClade10tips))
                bayestraitsdiscretedata[,1]<-names(horizontalClade10tips)
                bayestraitsdiscretedata[,2]<-as.numeric(as.character(horizontalClade10tips))
                bayestraitsdiscretedata[,2]<-as.numeric(bayestraitsdiscretedata[,2])-as.numeric(min(bayestraitsdiscretedata[,2]))
                bayestraitsdiscretedata[bayestraitsdiscretedata[,2]!=0,2]<-1
                bayestraitsdiscretedata[,3]<- as.numeric(as.character(Forestmember))
                bayestraitsdiscretedata<-as.data.frame(bayestraitsdiscretedata)
                
                command_vec3 <- c("2", "1") #option 1 = 2 discrete independent; option 2 = 1 maximum likelihood
                results_3 <- bayestraits(bayestraitsdiscretedata, currenttree, command_vec3)
                log_3 <- results_3$Log
                drift_results[counter,16]<-log_3$results$Lh
                
                command_vec4 <- c("3", "1") #option 1 = 3 discrete dependent; option 2 = 1 maximum likelihood
                results_4 <- bayestraits(bayestraitsdiscretedata, currenttree, command_vec4)
                log_4 <- results_4$Log
                drift_results[counter,17]<-log_4$results$Lh
                
              } #end of bayestraits option
              
              
              #For the simulations where there are more than two variants, we run the Bayestraits multistate option
              
              if(length(table(horizontalClade10tips))>2) {
                
                bayestraitsdiscretedata<-matrix(NA,nrow=length(horizontalClade10tips),ncol=2)
                colnames(bayestraitsdiscretedata)<-c("Species","Simulated")
                rownames(bayestraitsdiscretedata)<-c(1:length(horizontalClade10tips))
                bayestraitsdiscretedata[,1]<-names(horizontalClade10tips)
                bayestraitsdiscretedata[,2]<-as.numeric(as.character(horizontalClade10tips))
                bayestraitsdiscretedata[,2]<-as.numeric(bayestraitsdiscretedata[,2])-1
                bayestraitsdiscretedata<-as.data.frame(bayestraitsdiscretedata)
                
                
                command_vec1 <- c("1", "1") #option 1 = 1 multistate; option 2 = 1 maximum likelihood
                results_1 <- bayestraits(bayestraitsdiscretedata, currenttree, command_vec1)
                log_1 <- results_1$Log
                drift_results[counter,19]<-log_1$results$Lh
                
                command_vec2 <- c("1", "1", "Restrict q01 q02 q03 q10 q12 q13 q20 q21 q23 q30 q31 q32") #option 1 = 1 multistate; option 2 = 1 maximum likelihood; option restrict all dependent
                results_2 <- bayestraits(bayestraitsdiscretedata, currenttree, command_vec2)
                log_2 <- results_2$Log
                drift_results[counter,20]<-log_2$results$Lh
                
                
              } #end of bayestraits multistate option
              
            } #end of checking whether there are more than 1 variants in the sample
            
          } #end of checking whether there are more than 1 variants in Clade A
          
          counter<-counter+1
          
        } # end of the 10 random subsamples with 10% horizontal transmission within clade A
        
        
        
        # Horizontal transfer, 10% of societies change, adopt state most frequent within their own ecoregion
        
        for (randomsample in 1:1) {
          
          Forestsimulatedtips<-simulatedtips_neutral[names(simulatedtips_neutral) %in% names(Forestmember[Forestmember==1])]
          
          Foresttochange<-sample(Forestsimulatedtips[Forestsimulatedtips==names(table(Forestsimulatedtips)[2])],min(17,table(Forestsimulatedtips)[2]))
          levels(Foresttochange)<-c("1","2","3","4")
          Foresttochange[]<-1
          horizontalForest10tips<-simulatedtips_neutral
          horizontalForest10tips[names(horizontalForest10tips) %in% names(Foresttochange)]<-Foresttochange
          
          horizontalForest10tips<-droplevels(horizontalForest10tips)
          
          #Start filling in the data in the respective column
          drift_results[counter,1]<-tree_used
          ifelse(drift_variant<4,drift_results[counter,2]<-2,drift_results[counter,2]<-4)
          ifelse(drift_variant<4,drift_results[counter,3]<-drift_variant,drift_results[counter,3]<-drift_variant-3)
          drift_results[counter,4]<-repetition
          
          
          #Start with the full sample of 172 societies
          
          drift_results[counter,5]<-"horizontalForest10"
          drift_results[counter,6]<-length(unique(horizontalForest10tips))
          
          if (length(table(horizontalForest10tips))!=1) { 
            
            
            # the phylogenetic reconstruction with the drift model, the model under which the data were simulated
            equal_reconstruction <- ace(horizontalForest10tips, currenttree, type="discrete", model="ER")
            
            drift_results[counter,7]<-equal_reconstruction$loglik
            
            
            
            # the phylogenetic reconstruction with a one-state selection model
            transitionmatrix<-reconstruction_selection_four
            if(length(table(horizontalForest10tips))==3) transitionmatrix<-reconstruction_selection_three
            if(length(table(horizontalForest10tips))==2) transitionmatrix<-reconstruction_selection_two
            if(length(table(horizontalForest10tips))==1) transitionmatrix<-reconstruction_selection_one
            
            selection_to_a_reconstruction_twofold <- ace(horizontalForest10tips, currenttree, type="discrete", model=transitionmatrix)
            
            drift_results[counter,8]<-selection_to_a_reconstruction_twofold$loglik
            
            
            
            # the phylogenetic reconstruction with a pathway selection model
            transitionmatrix<-reconstruction_selection_four_pathway
            if(length(table(horizontalForest10tips))==3) transitionmatrix<-reconstruction_selection_three_pathway
            if(length(table(horizontalForest10tips))==2) transitionmatrix<-reconstruction_selection_two_pathway
            if(length(table(horizontalForest10tips))==1) transitionmatrix<-reconstruction_selection_one_pathway
            
            selection_to_a_reconstruction_pathway_slow <- fitDiscrete(dat=horizontalForest10tips, phy=currenttree, model="meristic")
            
            drift_results[counter,9]<-selection_to_a_reconstruction_pathway_slow$opt$lnL
            
            if(length(unique(horizontalForest10tips))>1) ifelse(names(table(horizontalForest10tips))[2]!=2, drift_results[counter,21]<-"TRUE",drift_results[counter,21]<-"FALSE") 
            
            
            # the phylogenetic reconstruction with an alternative pathway selection model
            transitionmatrix<-reconstruction_selection_four_pathway_toone
            if(length(table(horizontalForest10tips))==3) transitionmatrix<-reconstruction_selection_three_pathway_toone
            if(length(table(horizontalForest10tips))==2) transitionmatrix<-reconstruction_selection_two_pathway_toone
            if(length(table(horizontalForest10tips))==1) transitionmatrix<-reconstruction_selection_one_pathway_toone
            
            selection_to_a_reconstruction_pathway_toone <- ace(horizontalForest10tips, currenttree, type="discrete", model=transitionmatrix)
            
            drift_results[counter,10]<-selection_to_a_reconstruction_pathway_toone$loglik
            
            # the phylogenetic reconstruction with an alternative pathway selection model
            transitionmatrix<-reconstruction_selection_four_pathway_fromone
            if(length(table(horizontalForest10tips))==3) transitionmatrix<-reconstruction_selection_three_pathway_fromone
            if(length(table(horizontalForest10tips))==2) transitionmatrix<-reconstruction_selection_two_pathway_fromone
            if(length(table(horizontalForest10tips))==1) transitionmatrix<-reconstruction_selection_one_pathway_fromone
            
            selection_to_a_reconstruction_pathway_fromone <- ace(horizontalForest10tips, currenttree, type="discrete", model=transitionmatrix)
            
            drift_results[counter,18]<-selection_to_a_reconstruction_pathway_fromone$loglik
            
            
            
            # the phylogenetic reconstruction with an unconstrained model
            unconstrained_reconstruction <- ace(horizontalForest10tips, currenttree, type="discrete", model="ARD")
            
            drift_results[counter,11]<-unconstrained_reconstruction$loglik
            
            resultphylosiglambda<-phylosig(currenttree,as.numeric(horizontalForest10tips),method="lambda") 
            drift_results[counter,12]<-resultphylosiglambda$lambda
            resultphylosigK<-phylosig(currenttree,as.numeric(horizontalForest10tips),method="K") 
            drift_results[counter,13]<-resultphylosigK[1]
            
            
            #For the simulations where there are only two variants, we check whether by chance the variants
            #end up isolated in Clade A or associated with societies living in forest ecoregions
            #For this we use Bayestraits Discrete, assessing the likelihood of independent versus dependent evolution
            if(length(table(horizontalForest10tips))==2) {
              
              bayestraitsdiscretedata<-matrix(NA,nrow=length(horizontalForest10tips),ncol=3)
              colnames(bayestraitsdiscretedata)<-c("Species","Simulated","CladeA")
              rownames(bayestraitsdiscretedata)<-c(1:length(horizontalForest10tips))
              bayestraitsdiscretedata[,1]<-names(horizontalForest10tips)
              bayestraitsdiscretedata[,2]<-as.numeric(as.character(horizontalForest10tips))
              bayestraitsdiscretedata[,2]<-as.numeric(bayestraitsdiscretedata[,2])-as.numeric(min(bayestraitsdiscretedata[,2]))
              bayestraitsdiscretedata[bayestraitsdiscretedata[,2]!=0,2]<-1
              bayestraitsdiscretedata[,3]<- as.numeric(as.character(cladeAmember))
              bayestraitsdiscretedata<-as.data.frame(bayestraitsdiscretedata)
              
              command_vec1 <- c("2", "1") #option 1 = 2 discrete independent; option 2 = 1 maximum likelihood
              results_1 <- bayestraits(bayestraitsdiscretedata, currenttree, command_vec1)
              log_1 <- results_1$Log
              drift_results[counter,14]<-log_1$results$Lh
              
              command_vec2 <- c("3", "1") #option 1 = 3 discrete dependent; option 2 = 1 maximum likelihood
              results_2 <- bayestraits(bayestraitsdiscretedata, currenttree, command_vec2)
              log_2 <- results_2$Log
              drift_results[counter,15]<-log_2$results$Lh
              
              
              bayestraitsdiscretedata<-matrix(NA,nrow=length(horizontalForest10tips),ncol=3)
              colnames(bayestraitsdiscretedata)<-c("Species","Simulated","Forestmember")
              rownames(bayestraitsdiscretedata)<-c(1:length(horizontalForest10tips))
              bayestraitsdiscretedata[,1]<-names(horizontalForest10tips)
              bayestraitsdiscretedata[,2]<-as.numeric(as.character(horizontalForest10tips))
              bayestraitsdiscretedata[,2]<-as.numeric(bayestraitsdiscretedata[,2])-as.numeric(min(bayestraitsdiscretedata[,2]))
              bayestraitsdiscretedata[bayestraitsdiscretedata[,2]!=0,2]<-1
              bayestraitsdiscretedata[,3]<- as.numeric(as.character(Forestmember))
              bayestraitsdiscretedata<-as.data.frame(bayestraitsdiscretedata)
              
              command_vec3 <- c("2", "1") #option 1 = 2 discrete independent; option 2 = 1 maximum likelihood
              results_3 <- bayestraits(bayestraitsdiscretedata, currenttree, command_vec3)
              log_3 <- results_3$Log
              drift_results[counter,16]<-log_3$results$Lh
              
              command_vec4 <- c("3", "1") #option 1 = 3 discrete dependent; option 2 = 1 maximum likelihood
              results_4 <- bayestraits(bayestraitsdiscretedata, currenttree, command_vec4)
              log_4 <- results_4$Log
              drift_results[counter,17]<-log_4$results$Lh
              
            } #end of bayestraits option
            
            
            #For the simulations where there are more than two variants, we run the Bayestraits multistate option
            
            if(length(table(horizontalForest10tips))>2) {
              
              bayestraitsdiscretedata<-matrix(NA,nrow=length(horizontalForest10tips),ncol=2)
              colnames(bayestraitsdiscretedata)<-c("Species","Simulated")
              rownames(bayestraitsdiscretedata)<-c(1:length(horizontalForest10tips))
              bayestraitsdiscretedata[,1]<-names(horizontalForest10tips)
              bayestraitsdiscretedata[,2]<-as.numeric(as.character(horizontalForest10tips))
              bayestraitsdiscretedata[,2]<-as.numeric(bayestraitsdiscretedata[,2])-1
              bayestraitsdiscretedata<-as.data.frame(bayestraitsdiscretedata)
              
              
              command_vec1 <- c("1", "1") #option 1 = 1 multistate; option 2 = 1 maximum likelihood
              results_1 <- bayestraits(bayestraitsdiscretedata, currenttree, command_vec1)
              log_1 <- results_1$Log
              drift_results[counter,19]<-log_1$results$Lh
              
              command_vec2 <- c("1", "1", "Restrict q01 q02 q03 q10 q12 q13 q20 q21 q23 q30 q31 q32") #option 1 = 1 multistate; option 2 = 1 maximum likelihood; option restrict all dependent
              results_2 <- bayestraits(bayestraitsdiscretedata, currenttree, command_vec2)
              log_2 <- results_2$Log
              drift_results[counter,20]<-log_2$results$Lh
              
              
            } #end of bayestraits multistate option
            
          }
          
          counter<-counter+1
          
        } # end of the 10 random subsamples with 10% horizontal transmission within forests
        
        
        
        
        # Horizontal transfer, 10% of societies change, adopt state of a close neighbor
        
        for (randomsample in 1:1) {
          
          Latitudedistance<-as.matrix(dist(doubledata$lati))
          Longitudedistance<-as.matrix(dist(doubledata$long))
          Totaldistance<-Latitudedistance+Longitudedistance
          diag(Totaldistance)<-1000
          subsettochange<-sample(simulatedtips_neutral,24)
          
          for (findtheneighbor in 1:24) {
            
            findtherow<-Totaldistance[,colnames(Totaldistance)==names(subsettochange[findtheneighbor])]==min(Totaldistance[,colnames(Totaldistance)==names(subsettochange[findtheneighbor]) ] ,na.rm=T)
            
            
            subsettochange[findtheneighbor]<-simulatedtips_neutral[names(simulatedtips_neutral)==min(names(findtherow[findtherow==TRUE]))]
            
          }
          
          horizontalNeighbor10tips<-simulatedtips_neutral
          horizontalNeighbor10tips[names(horizontalNeighbor10tips) %in% names(subsettochange)]<-subsettochange
          
          horizontalNeighbor10tips<-droplevels(horizontalNeighbor10tips)
          
          #Start filling in the data in the respective column
          drift_results[counter,1]<-tree_used
          ifelse(drift_variant<4,drift_results[counter,2]<-2,drift_results[counter,2]<-4)
          ifelse(drift_variant<4,drift_results[counter,3]<-drift_variant,drift_results[counter,3]<-drift_variant-3)
          drift_results[counter,4]<-repetition
          
          
          #Start with the full sample of 172 societies
          
          drift_results[counter,5]<-"horizontalNeighbor10"
          drift_results[counter,6]<-length(unique(horizontalNeighbor10tips))
          
          if (length(table(horizontalNeighbor10tips))!=1) { 
            
            
            # the phylogenetic reconstruction with the drift model, the model under which the data were simulated
            equal_reconstruction <- ace(horizontalNeighbor10tips, currenttree, type="discrete", model="ER")
            
            drift_results[counter,7]<-equal_reconstruction$loglik
            
            
            
            # the phylogenetic reconstruction with a one-state selection model
            transitionmatrix<-reconstruction_selection_four
            if(length(table(horizontalNeighbor10tips))==3) transitionmatrix<-reconstruction_selection_three
            if(length(table(horizontalNeighbor10tips))==2) transitionmatrix<-reconstruction_selection_two
            if(length(table(horizontalNeighbor10tips))==1) transitionmatrix<-reconstruction_selection_one
            
            selection_to_a_reconstruction_twofold <- ace(horizontalNeighbor10tips, currenttree, type="discrete", model=transitionmatrix)
            
            drift_results[counter,8]<-selection_to_a_reconstruction_twofold$loglik
            
            
            
            # the phylogenetic reconstruction with a pathway selection model
            transitionmatrix<-reconstruction_selection_four_pathway
            if(length(table(horizontalNeighbor10tips))==3) transitionmatrix<-reconstruction_selection_three_pathway
            if(length(table(horizontalNeighbor10tips))==2) transitionmatrix<-reconstruction_selection_two_pathway
            if(length(table(horizontalNeighbor10tips))==1) transitionmatrix<-reconstruction_selection_one_pathway
            
            selection_to_a_reconstruction_pathway_slow <- fitDiscrete(dat=horizontalNeighbor10tips, phy=currenttree, model="meristic")
            
            drift_results[counter,9]<-selection_to_a_reconstruction_pathway_slow$opt$lnL
            
            if(length(unique(horizontalNeighbor10tips))>1) ifelse(names(table(horizontalNeighbor10tips))[2]!=2, drift_results[counter,21]<-"TRUE",drift_results[counter,21]<-"FALSE") 
            
            
            
            # the phylogenetic reconstruction with an alternative pathway selection model
            transitionmatrix<-reconstruction_selection_four_pathway_toone
            if(length(table(horizontalNeighbor10tips))==3) transitionmatrix<-reconstruction_selection_three_pathway_toone
            if(length(table(horizontalNeighbor10tips))==2) transitionmatrix<-reconstruction_selection_two_pathway_toone
            if(length(table(horizontalNeighbor10tips))==1) transitionmatrix<-reconstruction_selection_one_pathway_toone
            
            selection_to_a_reconstruction_pathway_toone <- ace(horizontalNeighbor10tips, currenttree, type="discrete", model=transitionmatrix)
            
            drift_results[counter,10]<-selection_to_a_reconstruction_pathway_toone$loglik
            
            # the phylogenetic reconstruction with an alternative pathway selection model
            transitionmatrix<-reconstruction_selection_four_pathway_fromone
            if(length(table(horizontalNeighbor10tips))==3) transitionmatrix<-reconstruction_selection_three_pathway_fromone
            if(length(table(horizontalNeighbor10tips))==2) transitionmatrix<-reconstruction_selection_two_pathway_fromone
            if(length(table(horizontalNeighbor10tips))==1) transitionmatrix<-reconstruction_selection_one_pathway_fromone
            
            selection_to_a_reconstruction_pathway_fromone <- ace(horizontalNeighbor10tips, currenttree, type="discrete", model=transitionmatrix)
            
            drift_results[counter,18]<-selection_to_a_reconstruction_pathway_fromone$loglik
            
            
            
            # the phylogenetic reconstruction with an unconstrained model
            unconstrained_reconstruction <- ace(horizontalNeighbor10tips, currenttree, type="discrete", model="ARD")
            
            drift_results[counter,11]<-unconstrained_reconstruction$loglik
            
            resultphylosiglambda<-phylosig(currenttree,as.numeric(horizontalNeighbor10tips),method="lambda") 
            drift_results[counter,12]<-resultphylosiglambda$lambda
            resultphylosigK<-phylosig(currenttree,as.numeric(horizontalNeighbor10tips),method="K") 
            drift_results[counter,13]<-resultphylosigK[1]
            
            
            #For the simulations where there are only two variants, we check whether by chance the variants
            #end up isolated in Clade A or associated with societies living in forest ecoregions
            #For this we use Bayestraits Discrete, assessing the likelihood of independent versus dependent evolution
            if(length(table(horizontalNeighbor10tips))==2) {
              
              bayestraitsdiscretedata<-matrix(NA,nrow=length(horizontalNeighbor10tips),ncol=3)
              colnames(bayestraitsdiscretedata)<-c("Species","Simulated","CladeA")
              rownames(bayestraitsdiscretedata)<-c(1:length(horizontalNeighbor10tips))
              bayestraitsdiscretedata[,1]<-names(horizontalNeighbor10tips)
              bayestraitsdiscretedata[,2]<-as.numeric(as.character(horizontalNeighbor10tips))
              bayestraitsdiscretedata[,2]<-as.numeric(bayestraitsdiscretedata[,2])-as.numeric(min(bayestraitsdiscretedata[,2]))
              bayestraitsdiscretedata[bayestraitsdiscretedata[,2]!=0,2]<-1
              bayestraitsdiscretedata[,3]<- as.numeric(as.character(cladeAmember))
              bayestraitsdiscretedata<-as.data.frame(bayestraitsdiscretedata)
              
              command_vec1 <- c("2", "1") #option 1 = 2 discrete independent; option 2 = 1 maximum likelihood
              results_1 <- bayestraits(bayestraitsdiscretedata, currenttree, command_vec1)
              log_1 <- results_1$Log
              drift_results[counter,14]<-log_1$results$Lh
              
              command_vec2 <- c("3", "1") #option 1 = 3 discrete dependent; option 2 = 1 maximum likelihood
              results_2 <- bayestraits(bayestraitsdiscretedata, currenttree, command_vec2)
              log_2 <- results_2$Log
              drift_results[counter,15]<-log_2$results$Lh
              
              
              bayestraitsdiscretedata<-matrix(NA,nrow=length(horizontalNeighbor10tips),ncol=3)
              colnames(bayestraitsdiscretedata)<-c("Species","Simulated","Forestmember")
              rownames(bayestraitsdiscretedata)<-c(1:length(horizontalNeighbor10tips))
              bayestraitsdiscretedata[,1]<-names(horizontalNeighbor10tips)
              bayestraitsdiscretedata[,2]<-as.numeric(as.character(horizontalNeighbor10tips))
              bayestraitsdiscretedata[,2]<-as.numeric(bayestraitsdiscretedata[,2])-as.numeric(min(bayestraitsdiscretedata[,2]))
              bayestraitsdiscretedata[bayestraitsdiscretedata[,2]!=0,2]<-1
              bayestraitsdiscretedata[,3]<- as.numeric(as.character(Forestmember))
              bayestraitsdiscretedata<-as.data.frame(bayestraitsdiscretedata)
              
              command_vec3 <- c("2", "1") #option 1 = 2 discrete independent; option 2 = 1 maximum likelihood
              results_3 <- bayestraits(bayestraitsdiscretedata, currenttree, command_vec3)
              log_3 <- results_3$Log
              drift_results[counter,16]<-log_3$results$Lh
              
              command_vec4 <- c("3", "1") #option 1 = 3 discrete dependent; option 2 = 1 maximum likelihood
              results_4 <- bayestraits(bayestraitsdiscretedata, currenttree, command_vec4)
              log_4 <- results_4$Log
              drift_results[counter,17]<-log_4$results$Lh
              
            } #end of bayestraits option
            
            
            
            #For the simulations where there are more than two variants, we run the Bayestraits multistate option
            
            if(length(table(horizontalNeighbor10tips))>2) {
              
              bayestraitsdiscretedata<-matrix(NA,nrow=length(horizontalNeighbor10tips),ncol=2)
              colnames(bayestraitsdiscretedata)<-c("Species","Simulated")
              rownames(bayestraitsdiscretedata)<-c(1:length(horizontalNeighbor10tips))
              bayestraitsdiscretedata[,1]<-names(horizontalNeighbor10tips)
              bayestraitsdiscretedata[,2]<-as.numeric(as.character(horizontalNeighbor10tips))
              bayestraitsdiscretedata[,2]<-as.numeric(bayestraitsdiscretedata[,2])-1
              bayestraitsdiscretedata<-as.data.frame(bayestraitsdiscretedata)
              
              
              command_vec1 <- c("1", "1") #option 1 = 1 multistate; option 2 = 1 maximum likelihood
              results_1 <- bayestraits(bayestraitsdiscretedata, currenttree, command_vec1)
              log_1 <- results_1$Log
              drift_results[counter,19]<-log_1$results$Lh
              
              command_vec2 <- c("1", "1", "Restrict q01 q02 q03 q10 q12 q13 q20 q21 q23 q30 q31 q32") #option 1 = 1 multistate; option 2 = 1 maximum likelihood; option restrict all dependent
              results_2 <- bayestraits(bayestraitsdiscretedata, currenttree, command_vec2)
              log_2 <- results_2$Log
              drift_results[counter,20]<-log_2$results$Lh
              
              
            } #end of bayestraits multistate option
            
          }
          
          counter<-counter+1
          
          
        } # end of the 10 random subsamples with 10% horizontal transmission from a close neighbor
        
        print(counter)
        print(repetition)
        print(tree_used)
        print(drift_variant)
        
        
      } #only need to run if there is more than one variant  
      
    } # end of the 1000 repetition loop
    
    
  }  # end of the drift variants model
  
  
} #end of the tree variant loop

options(warn = oldw)

#Store the output file on your computer
write.csv(drift_results[1:counter,],file="SimulateCulture_DriftResults.csv")
