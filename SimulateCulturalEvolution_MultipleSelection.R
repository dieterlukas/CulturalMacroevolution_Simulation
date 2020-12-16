

#------------------------------------------------------------------------------------------
# Simulations for the manuscript - selection
#------------------------------------------------------------------------------------------

# These simulations and their outcomes are described in:
# Lukas, D., Towner, M., & Mulder, M. B. (2020). The Potential to Infer the Historical Pattern of Cultural Macroevolution as Illustrated by the Western North American Indian Societies.
# https://osf.io/preprints/socarxiv/tjvgy/



#------------------------------------------------------------------------------------------
# This set of simulations and analyses aims to assess the potential false negative error rate 
# of phylogenetic reconstructions of the historical evolution of cultural traits

# The code simulates changes in an arbitrary cultural trait across a known phylogeny, assuming that  
# the rate of change differs in different parts of the tree. This reflects a trait that experienced
# different selective regimes, reflecting either lineage-specific evolution or selection associated
# with differences in ecology.
# The resulting distribution of the different states of the trait are analysed with a variety of methods
# that are frequently used in biological and cultural phylogenetics. A false negative error occurs
# if such an analysis wrongly infers that one of the variants of the trait has not been under selection.


# In addition, the simulation assesses how frequently an analyses would wrongly infer that the distribution
# of states of the trait differs between societies in different lineages or with different ecologies, even
# though the changes are simulated to occur equally throughout the tree.

# The analyses examine three potential errors:
# 1) wrongly identifying that selection did not occur on the trait:
#         determines whether an analysis wrongly supports the inference that all changes between the states
#         of the trait have occurred at the same rate - even though the simulation is set up
#         such that changes to and from some states are more likely to occur.
# 2) wrongly missing lineage differences in changes in the trait:
#         determines whether an analysis wrongly supports the inference that changes between states
#         occurred equally likely in different parts of the tree, even though the simulation is set up
#         such that changes occur at different rates in some lineages or in some ecologies


# The setup includes further modifications to assess whether the false negative error rates are influenced by
# i)   properties of the sample:
#         assesses how the inferences change if information is available for only a subset of tree (subclade),
#         or if information is missing from societies throughout the tree
# ii)  properties of the tree:
#         assesses how the inferences change if known branch lengths are used, if branch lengths are modified
#         to assume that changes in the trait are associated with splits between societies but not necessarily 
#         happen constantly over time, and if phylogeny is imbalanced towards early or late diversification
# iii) properties of the trait:
#         assesses how inferences change if the trait is split into more or fewer states, 
#         if changes occur frequently or infrequently, and if states can be transmitted horizontally


# The simulations focus on traits that can be split into two or more discrete categories.
# This reflects a large number of cultural traits, which are either classified as present/absent or
# into different separate categories. The simulations currently do not reflect continuous traits.

# A matching set of simulations and analyses to assess the potential false positive rate
# can be found here: https://github.com/dieterlukas/CulturalMacroevolution_Simulation/blob/master/SimulateCulturalEvolution_Drift.R

# This code can be adapted to assess the potential risk of false negative errors in a user-provided dataset.
# For this, you need a phylogeny in nexus format, and, if relevant, information on an ecological variable
# for the societies in the sample. Load these instead of the example files and adjust the steps after line 130.

#------------------------------------------------------------------------------------------



#Load necessary libraries
# The code relies on a number of phylogenetic packages to manipulate and analyse the data
library(ape)
library(geiger)
library(phytools)
library(OUwie)

# The code also uses some packages that facilitate loading and structuring the data
library(dplyr)
library(readr)

# Finally, the simulation relies on the software Bayestraits (http://www.evolution.rdg.ac.uk/SoftwareMain.html)
# It starts the software from within R using the package Bayestraits Wrapper
# To install the package, follow the instructions by the developer (http://www.randigriffin.com/projects/btw.html)
# It is also necessary that you have Bayestraits on your computer, and that you set the working directory in R to a location with a copy of Bayestraits
library(btw)



#------------------------------------------------------------------------------------------
# There are two sets of societies here as example:
# Western North American Indigeneous societies (WNAI) and societies from Australia who speak a language of the Pama Nyungan family

# To run the simulations across a phylogeny of Western North American Indigeneous societies, run this next line
Option<-"WNAI"

# To run the simulations across a phylogeny of Pama PamaNyungan speaking societies, run this next line
Option<-"PamaNyungan"


#------------------------------------------------------------------------------------------
# Load relevant phylogenies from GitHub


#Load WNAI tree constructed from the language classifications
if(Option=="WNAI") {
  Americantree<-read.nexus(url("https://raw.githubusercontent.com/dieterlukas/CulturalMacroevolution_Simulation/master/WNAI_tree_forsimulation.nex"))
}

if(Option=="PamaNyungan") {
  #Load Pama Nyungan tree provided by Bouckaert et al. 2018 https://doi.org/10.1038/s41559-018-0489-3
  PamaNyungantree<-read.nexus(url("https://raw.githubusercontent.com/dieterlukas/CulturalMacroevolution_Simulation/master/PamaNyungan_tree_forsimulation.nex"))
}


#------------------------------------------------------------------------------------------
# Load associated data from GitHub

# For both sets of societies, there is information on one ecological variable. 

# For the WNAI, the variable classifies the ecoregion that each society predominantly uses (32 different ecosystems; the data are from a version of the WNAI dataset provided by Dow & Eff http://intersci.ss.uci.edu/wiki/index.php/Materials_for_cross-cultural_research)
if(Option=="WNAI") {
  WNAIdata<-read_csv(url("https://raw.githubusercontent.com/dieterlukas/CulturalMacroevolution_Simulation/master/WNAI_data_forsimulation.csv"))
  WNAIdata<-data.frame(WNAIdata)
  
  # In addition, we load the geographic location for each society (latitude/longitude) for the horizontal transmission among neighbors; the data are from the same dataset
  WNAIlocations<-read_csv(url("https://raw.githubusercontent.com/dieterlukas/CulturalMacroevolution_Simulation/master/WNAI_locations_forsimulation.csv"))
  WNAIlocations<-data.frame(WNAIlocations)
  
}

# For the Pama Nyungan, the variable classifies the main mode of subsistence (hunter-gatherer versus food produce; the data are from Derungs et al. 2018 https://github.com/curdon/linguisticDensity_ProcB_derungsEtAl)
if(Option=="PamaNyungan") {
  PamaNyungandata<-read_csv(url("https://raw.githubusercontent.com/dieterlukas/CulturalMacroevolution_Simulation/master/PamaNyungan_data_forsimulation.csv"))
  PamaNyungandata<-data.frame(PamaNyungandata)
  rownames(PamaNyungandata)<-PamaNyungandata$society
  
  # In addition, we load the geographic location for each society (latitude/longitude) for the horizontal transmission among neighbors; the data are from the publication that describes the phylogeny, based on the Bowern, Claire (2016). Chirila: Contemporary and Historical Resources for the Indigenous Languages of Australia. Language Documentation and Conservation. Vol 10. http://www.pamanyungan.net/chirila/
  PamaNyunganlocations<-read_csv(url("https://raw.githubusercontent.com/dieterlukas/CulturalMacroevolution_Simulation/master/PamaNyungan_locations_forsimulation.csv"))
  PamaNyunganlocations<-data.frame(PamaNyunganlocations)
  rownames(PamaNyunganlocations)<-PamaNyunganlocations$Language
}








#------------------------------------------------------------------------------------------

#------------------------------------------------------------------------------------------
# One part of the simulation is to assess whether the shape of the tree, and in particular the branch lenghts have an influence on the inferences
# For this, we built four additional variants of each phylogenetic tree:
# Grafentree: a tree with branch lengths based on Grafen's method (all tips equidistant from root, branch length depends onnumber of nodes between root and tip)
# Onetree: a tree with all branch lengths set to have the same length of one
# Earlytree: a tree with early diversification and long branches leading to the tips
# Latetree: a tree with recent diversification and long branches between clades


if(Option=="WNAI") {
  #Modify the tree of the WNAI societies for the analyses
  
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
}


#------------------------------------------------------------------------------------------

if(Option=="PamaNyungan") {
  #Modify the tree of the PamaNyungan societies for the analyses
  
  #Add branch lengths to the tree based on Grafen's method (all tips equidistant from root, branch length depends onnumber of nodes between root and tip)
  Grafentree<-compute.brlen(PamaNyungantree,method="Grafen")
  
  #Add branch lengths to the tree assuming that all branches have the same length of one
  Onetree<-compute.brlen(PamaNyungantree,1)
  
  #Add branch lengths to the tree with early diversification and long branches to the tips
  Earlytree<-compute.brlen(PamaNyungantree,method="Grafen",power=0.25)
  
  #Add branch lengths to the tree with lots of recent diversification and long branches between clades
  Latetree<-compute.brlen(PamaNyungantree,method="Grafen",power=1.5)
  
  #Some analyses need a rooted, fully bifurcating tree
  PamaNyungantree<-root(PamaNyungantree,node=307)
  PamaNyungantree<-multi2di(PamaNyungantree)
  
  Grafentree<-root(Grafentree,node=307)
  Grafentree<-multi2di(Grafentree)
  Grafentree<-compute.brlen(Grafentree,method="Grafen")
  
  Onetree<-root(Onetree,node=307)
  Onetree<-multi2di(Onetree)
  Onetree<-compute.brlen(Onetree,1)
  
  Latetree<-root(Latetree,node=307)
  Latetree<-multi2di(Latetree)
  Latetree<-compute.brlen(Latetree,method="Grafen",power=1.5)
  
  Earlytree<-root(Earlytree,node=307)
  Earlytree<-multi2di(Earlytree)
  Earlytree<-compute.brlen(Earlytree,method="Grafen",power=0.25)
  
}


#------------------------------------------------------------------------------------------

#------------------------------------------------------------------------------------------
# Another part of the simulations is to assess whether the size of the phylogeny (number of societies) has an influence on inferences
# For this, we create subsets from each tree reflecting major clades to get smaller trees with more limited numbers of societies

# An additional part of the simulations is to determine whether changes occurred differently in different parts of the tree
# Also for this, we want to classify the tree and societies into different clades that reflect major lineages

#Identify the major subclades in the tree


if(Option=="WNAI") {
  # For the WNAI phylogeny
  
  # We identify the nodes that are at the base of each of the clades
  # For the WNAI, clades are based on higher-order language classifications
  # We pick two societies with the same classification from opposite ends of the lineage to identify the ancestral node
  cladeA_basalnode<-mrca(Americantree)["92","110"]  #basal node 192
  cladeB_basalnode<-mrca(Americantree)["3","1"] #basal node 174
  cladeC_basalnode<-mrca(Americantree)["172","72"] #basal node 240
  cladeD_basalnode<-mrca(Americantree)["6","160"] #basal node 193
  cladeE_basalnode<-mrca(Americantree)["44","119"] #basal node 227
  cladeF_basalnode<-mrca(Americantree)["162","10"] #basal node 211
  
  #Create subsets of trees according to the major language groups based on the previously identified basal nodes
  #We create lists with the names of the societies in each clade (each society that is descended from the basal node)
  cladeA <- tips(Americantree, cladeA_basalnode) # 102 societies "NorthernAmerind"
  cladeB <- tips(Americantree, cladeB_basalnode) # 28 societies "NaDene"
  cladeC <- tips(Americantree, cladeC_basalnode) # 42 societies "CentralAmerind"
  cladeD <- tips(Americantree, cladeD_basalnode) # 36 societies "Penutian"
  cladeE <- tips(Americantree, cladeE_basalnode) # 27 societies "Hokan"
  cladeF <- tips(Americantree, cladeF_basalnode) # 39 societies "Almosan"
  
  #We create smaller trees representing each of these clades
  cladeAtree<-keep.tip(Americantree,cladeA)
  cladeBtree<-keep.tip(Americantree,cladeB)
  cladeCtree<-keep.tip(Americantree,cladeC)
  cladeDtree<-keep.tip(Americantree,cladeD)
  cladeEtree<-keep.tip(Americantree,cladeE)
  cladeFtree<-keep.tip(Americantree,cladeF)
}


if(Option=="PamaNyungan") {
  # For the PamaNyungan phylogeny
  
  # We identify the nodes that are at the base of each of the clades
  # We pick two societies from within the same clade from opposite ends of the lineage to identify the ancestral node
  cladeA_basalnode<-mrca(PamaNyungantree)["WangkumaraMcDWur","Adnyamathanha"]  #309
  cladeB_basalnode<-mrca(PamaNyungantree)["Minkin","Zorc"] #basal node 508
  cladeC_basalnode<-mrca(PamaNyungantree)["Kunjen","Wulguru"] #basal node 424
  cladeD_basalnode<-mrca(PamaNyungantree)["Piangil","Iyora"] #basal node 528
  cladeE_basalnode<-mrca(PamaNyungantree)["KKY","Olkola"] #basal node 425
  cladeF_basalnode<-mrca(PamaNyungantree)["Nyamal","Wajarri"] #basal node 324
  
  
  #Create subsets of trees according to the major language groups based on the previously identified basal nodes
  #We create lists with the names of the societies in each clade (each society that is descended from the basal node)
  cladeA <- tips(PamaNyungantree, cladeA_basalnode) # 113 societies "WesternExpansion"
  cladeB <- tips(PamaNyungantree, cladeB_basalnode) # 21 societies "NorthernCape"
  cladeC <- tips(PamaNyungantree, cladeC_basalnode) # 78 societies "PamaMaric"
  cladeD <- tips(PamaNyungantree, cladeD_basalnode) # 81 societies "SouthernCape"
  cladeE <- tips(PamaNyungantree, cladeE_basalnode) # 48 societies "NorthEast"
  cladeF <- tips(PamaNyungantree, cladeF_basalnode) # 34 societies "Southwest"
  
  #We create smaller trees representing each of these clades
  cladeAtree<-keep.tip(PamaNyungantree,cladeA)
  cladeBtree<-keep.tip(PamaNyungantree,cladeB)
  cladeCtree<-keep.tip(PamaNyungantree,cladeC)
  cladeDtree<-keep.tip(PamaNyungantree,cladeD)
  cladeEtree<-keep.tip(PamaNyungantree,cladeE)
  cladeFtree<-keep.tip(PamaNyungantree,cladeF)
  
}



#------------------------------------------------------------------------------------------

#------------------------------------------------------------------------------------------
# For the analyses to determine whether changes occurred differently in different parts of the tree, there are two specifications
# The first one takes the information from above to contrast changes in the largest clade (CladeA) from changes in the remaining part of the tree
# The second one takes information on an ecological variables to classify societies into two groups: 
# for the American societies this is whether they primarily use a forest habitat or not
# for the PamaNyungan societies this is whether they are hunter-gatherers or food-producers

# To determine whether changes occurred differently, we need to match this information from the societies to the tree
# For the clade-based distinction, this labels all branches within the clade one way an all other branches another way
# For the ecology-based distinction, branches are labelled based on a phylogenetic reconstruction of the most likely history of the ecological variable


if(Option=="WNAI") {
  #For the analyses asking whether selection operated differently in different clades, 
  #we split the tree according to the largest clades - CladeA versus the remaining societies
  
  #Create dummy variables designating subclade membership for each society 
  WNAIdata$cladeA<-rep(0,172)
  WNAIdata[WNAIdata$tribes %in% cladeA,]$cladeA<-1
  cladeAmember<-WNAIdata$cladeA
  names(cladeAmember)<-WNAIdata$tribes
  
  #Label all the nodes of the respective subclades
  simmapattemptCladeA<-make.simmap(Grafentree,cladeAmember)
  
  #We also want that selection operates differently in one of the ecoregions - all of them that are forests
  # Information on the ecoregion is available for all societies in the sample. 
  WNAIdata$Forests<-as.numeric(grepl("Forest",WNAIdata$ecoregion))
  Forestmember<-WNAIdata$Forests
  names(Forestmember)<-WNAIdata$tribes
  Ecologymember<-Forestmember
  
  #Label all the nodes according to whether the ancestral society was likely to be associated with forests
  #This uses the stochastic mapping approach for binary characters as implemented in phytools
  simmapattemptEcology<-make.simmap(Grafentree,Ecologymember)
  
}

if(Option=="PamaNyungan") {
  #For the analyses asking whether selection operated differently in different clades, 
  #we split the tree according to the largest clades - CladeA versus the remaining societies
  
  #Create dummy variables designating subclade membership for each society 
  PamaNyungandata$cladeA<-rep(0,306)
  PamaNyungandata[PamaNyungandata$society %in% cladeA,]$cladeA<-1
  cladeAmember<-PamaNyungandata$cladeA
  names(cladeAmember)<-PamaNyungandata$society
  
  #Label all the nodes of the respective subclades
  simmapattemptCladeA<-make.simmap(Grafentree,cladeAmember)
  
  #We also want that selection operates differently in one of the ecologies - all of them that are hunter gatherers
  PamaNyungandata$huntergatherer<-as.numeric(PamaNyungandata$subsistence=="huntergatherer")
  Huntergatherermember<-PamaNyungandata$huntergatherer
  names(Huntergatherermember)<-PamaNyungandata$society
  Ecologymember<-Huntergatherermember
  
  imputedecology<-phylo.impute(PamaNyungantree,as.matrix(Ecologymember))
  
  Ecologymember<-as.numeric(imputedecology[,1]>0.5)
  names(Ecologymember)<-names(imputedecology[,1])
  
  #We do not have information on the subsistence mode for all societies in the sample. 
  # For those with missing data, we assign an equal probability that they are either hunter-gatherers or food producers
  Huntergatherermembermatrix<-to.matrix(Huntergatherermember,c("0","1"))
  Huntergatherermembermatrix[is.na(Huntergatherermember),]<-c(0.5,0.5)
  
  #Label all the nodes according to whether ancestral society was likely to be hunting/gathering
  #This uses the stochastic mapping approach for binary characters as implemented in phytools
  simmapattemptEcology<-make.simmap(Grafentree,Huntergatherermembermatrix)
  
  #paintedPamaNyungantree<-paintBranches(PamaNyungantree,edge=PamaNyungantree$edge[simmapattemptEcology$mapped.edge[,1]>0.5,2],"1","0")
  
}





# This completes the preparation of the tree and data
#------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------




#------------------------------------------------------------------------------------------

#------------------------------------------------------------------------------------------
# The following part sets up matrices reflecting the assumed transitions for the simulated traits
# This simulation assumes that changes in the simulated trait occur with a specified rate over time
# The different variants make some changes more likely and others less likely (or impossible)
# It assumes that transitions occur the same everywhere across the tree
# There are nine matrices: assuming either two or four states of the trait; 
# slow, medium, or fast rates of transition; and either favoring one trait or favoring a pathway of change.

#Create the matrices for the simulations of the data on the phylogenetic tree

#Changes to b only occur in Ecology 1 and they occur very fast
Q_onlyandfastInEcology1<-setNames(list(matrix(c(-2,0,12,-2),2,2,dimnames=list(letters[1:2],letters[1:2])),matrix(c(0,0,0,0),2,2,dimnames=list(letters[1:2],letters[1:2]))),c("0","1"))

#Changes to b only occur in Ecology 1 and they occur relatively slowly
Q_onlyandslowInEcology1<-setNames(list(matrix(c(-1,0,2,-1),2,2,dimnames=list(letters[1:2],letters[1:2])),matrix(c(0,0,0,0),2,2,dimnames=list(letters[1:2],letters[1:2]))),c("0","1"))

#Changes to b occur relatively slow and are twice as likely in Ecology 1
Q_moreandslowInEcology1<-setNames(list(matrix(c(-1,1,4,-1),2,2,dimnames=list(letters[1:2],letters[1:2])),matrix(c(-2,2,2,-2),2,2,dimnames=list(letters[1:2],letters[1:2]))),c("0","1"))

#Changes to b occur relatively slow and are four times as likely in Ecology 1
Q_moreandfastInEcology1<-setNames(list(matrix(c(-1,0.5,8,-1),2,2,dimnames=list(letters[1:2],letters[1:2])),matrix(c(-2,2,2,-2),2,2,dimnames=list(letters[1:2],letters[1:2]))),c("0","1"))







#-------------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------




# The previous sections were to prepare the data and set up the simulations
# We now select the correct settings for the specific run of the simulations

# The following selects the dataset, and sets up loops to repeat simulations 
# across all variants of the trees and using the different drift models


if(Option=="WNAI") {data<-WNAIdata}
if(Option=="PamaNyungan") {data<-PamaNyungandata}
if(Option=="WNAI") {locationdata<-WNAIlocations}
if(Option=="PamaNyungan") {locationdata<-PamaNyunganlocations}

if(Option=="WNAI") {tree_variants <- c("Grafentree","Onetree","Earlytree","Latetree")}
if(Option=="PamaNyungan") {tree_variants <- c("PamaNyungantree","Grafentree","Onetree","Earlytree","Latetree")}


selection_variants<-c("Q_onlyandfastInEcology1","Q_onlyandslowInEcology1","Q_moreandslowInEcology1","Q_moreandfastInEcology1")

counter<-1

# To make it easier to follow the process of the simulation we suppress warnings
oldw <- getOption("warn")
options(warn = -1)

#Change the number of repetitions to modify the number of indepedent simulations with a given drift model on a given tree are being performed and analysed
#There are many different independent analyses for each simulation, so the process takes a significant amount of time for large numbers of repetitions of simulations
repetitions<-10

#We prepare a data frame to store all the output
#Each row will have the results from a single simulation (with a given simulation model on a given tree)
#The columns contain basic information on the settings and the output from the various reconstruction methods
selection_results<-matrix(data=NA,nrow=500000,ncol=22)
colnames(selection_results)<-c("Tree","NumberOfVariantsInModel","RateOfChange","Repetition","Sample","NumberOfVariantsObserved","loglik_drift","loglik_selectiononestate","loglik_pathways","loglik_pathwaystoone","loglik_unconstrained","phylosigLambda","phylosigK","BayestraitsIndependent_CladeA","BayestraitsDependent_CladeA","BayestraitsIndependent_Forest","BayestraitsDependent_Forest","loglik_pathwaysfromone","BayestraitsMultistate_Free","BayestraitsMultistate_Equal","Pathway_impossible","selection_variant")






#-------------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------

# The simulation loops over the four (WNAI) or five (PamaNyungan) different phylogenies; for each phylogeny
# it simulates trait histories using one of the six drift variants; meaning there are 24/30 settings. 
# The number of repetitions set above apply to each setting (so with 10 repetitions there are 240/300 loops)

# For each sample, setting and repetition, the following steps occur:
# - first the evolutionary history of a trait is modelled across the phylogeny
# - next, different samples are generated based on this evolutionary history
#       the full sample of all societies; 
#       six samples reflecting clades A-F;  
#       a sample where data for 25% of randomly drawn societies is missing
#       a sample where data for 50% of randomly drawn societies is missing
#       a sample where data is missing for 50% of societies with the least frequent state
#       a sample where data is missing for 50% of societies with the most frequent state
#       a sample where data is missing for 50% of societies with the dominant ecology
#       a sample where data is missing for 50% of societies with the most frequent state
#       a sample where horizontal transfer among contemporary societies happened, with 10% of societies adopting a random other state
#       a sample where horizontal transfer among contemporary societies happened, with 10% of societies adopting the state that is most frequent in their clade
#       a sample where horizontal transfer among contemporary societies happened, with 10% of societies adopting the state of their close geographic neighbor


# For each sample (within each setting and repetition), the following analyses are performed:
# - next, the simulated distribution of states of the trait across societes is used as input data in 
#   phylogenetic reconstructions with drift model, one of the four selection models, and an unconstrained model
# - next, the simulated distribution of states of the trait across societies is passed on to Bayestraits;
#   when there are only two states of the trait, it investigates whether the distribution of the traits is
#   different between clade A and the remainder of the tree and different between societies with the two types
#   of ecology; when there are three or more states it determines whether all transitions are equally likely or not.








#This will start the loop; 
#because of the large number of analyses for each simulation (different techniques, various subsamples), this takes noticeable computer time

for (condition_variant in 1:2) {

for (tree_variant in 1:length(tree_variants)) {
  
  tree_used<-tree_variants[tree_variant]
  if(tree_used=="PamaNyungantree") currenttree<-PamaNyungantree
  if(tree_used=="Grafentree") currenttree<-Grafentree
  if(tree_used=="Onetree") currenttree<-Onetree
  if(tree_used=="Earlytree") currenttree<-Earlytree
  if(tree_used=="Latetree") currenttree<-Latetree
  
  for (selection_variant in 1:length(selection_variants)) {
    
    selectionvariant_used<-selection_variants[selection_variant]
    if(selectionvariant_used=="Q_onlyandfastInEcology1") current_selectionvariant<-Q_onlyandfastInEcology1
    if(selectionvariant_used=="Q_moreandslowInEcology1") current_selectionvariant<-Q_moreandslowInEcology1
    if(selectionvariant_used=="Q_onlyandslowInEcology1") current_selectionvariant<-Q_onlyandslowInEcology1
    if(selectionvariant_used=="Q_moreandfastInEcology1") current_selectionvariant<-Q_moreandfastInEcology1
        

    for (repetition in 1:repetitions) {
      
      
      #Simulate the data with the specified tree and the specified drift model
      
      if(condition_variant==1) {simulatedtips_neutral<-sim.multiMk(tree=simmapattemptCladeA,Q=current_selectionvariant,anc="a")}
      if(condition_variant==2) {simulatedtips_neutral<-sim.multiMk(tree=simmapattemptEcology,Q=current_selectionvariant,anc="a")}
      
      
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
        selection_results[counter,1]<-tree_used
        selection_results[counter,2]<-4
        if(selection_variant<4) selection_results[counter,2]<-2
        selection_results[counter,3]<-selection_variant-3
        if(selection_variant<4) selection_results[counter,3]<-selection_variant
        selection_results[counter,4]<-repetition
        selection_results[counter,22]<-selectionvariant_used
        
        
        #Start the phylogenetic reconstruction with the simulated data
        
        
        #Start with the full sample of societies
        
        selection_results[counter,5]<-"full"
        selection_results[counter,6]<-length(unique(simulatedtips_neutral))
        
        
        numeric_simulatedtipds_neutral<-as.numeric(simulatedtips_neutral)
        names(numeric_simulatedtipds_neutral)<-names(simulatedtips_neutral)
        resultphylosiglambda<-phylosig(currenttree,numeric_simulatedtipds_neutral,method="lambda") 
        selection_results[counter,12]<-resultphylosiglambda$lambda
        resultphylosigK<-phylosig(currenttree,numeric_simulatedtipds_neutral,method="K") 
        selection_results[counter,13]<-resultphylosigK[1]
        
        
        #For the simulations where there are only two variants, we check whether by chance the variants
        #end up isolated in Clade A or associated with societies having different ecologies
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
          selection_results[counter,14]<-log_1$results$Lh
          
          command_vec2 <- c("3", "1") #option 1 = 3 discrete dependent; option 2 = 1 maximum likelihood
          results_2 <- bayestraits(bayestraitsdiscretedata, currenttree, command_vec2)
          log_2 <- results_2$Log
          selection_results[counter,15]<-log_2$results$Lh
          
          
          bayestraitsdiscretedata<-matrix(NA,nrow=length(simulatedtips_neutral),ncol=3)
          colnames(bayestraitsdiscretedata)<-c("Species","Simulated","Ecologymember")
          rownames(bayestraitsdiscretedata)<-c(1:length(simulatedtips_neutral))
          bayestraitsdiscretedata[,1]<-names(simulatedtips_neutral)
          bayestraitsdiscretedata[,2]<-as.numeric(as.character(simulatedtips_neutral))
          bayestraitsdiscretedata[,2]<-as.numeric(bayestraitsdiscretedata[,2])-as.numeric(min(bayestraitsdiscretedata[,2]))
          bayestraitsdiscretedata[bayestraitsdiscretedata[,2]!=0,2]<-1
          bayestraitsdiscretedata[,3]<- as.numeric(as.character(Ecologymember))
          bayestraitsdiscretedata[,3]<-as.numeric(as.character(bayestraitsdiscretedata[,3]))
          bayestraitsdiscretedata<-as.data.frame(bayestraitsdiscretedata)
          bayestraitsdiscretedata[is.na(bayestraitsdiscretedata[,3]),3]<-"-"
          command_vec3 <- c("2", "1") #option 1 = 2 discrete independent; option 2 = 1 maximum likelihood
          results_3 <- bayestraits(bayestraitsdiscretedata, currenttree, command_vec3)
          log_3 <- results_3$Log
          selection_results[counter,16]<-log_3$results$Lh
          
          command_vec4 <- c("3", "1") #option 1 = 3 discrete dependent; option 2 = 1 maximum likelihood
          results_4 <- bayestraits(bayestraitsdiscretedata, currenttree, command_vec4)
          log_4 <- results_4$Log
          selection_results[counter,17]<-log_4$results$Lh
          
        } #end of bayestraits 2 option
        
      
        
        
        
        
        
        
        counter<-counter+1
        
        
        
      
        
        
        #--#--#--#--#--#--#        #-#-#-#-#-#-#-#-#-#-#-#
        
        # Subset the tree and the data to assume that 25% of data is randomly missing
        
        
        
        ThreeQuarterClade <- currenttree$tip.label
        ThreeQuarterClade <- sample(ThreeQuarterClade)
        ThreeQuarterClade <- ThreeQuarterClade[1:129]
        
        ThreeQuartertree<-keep.tip(currenttree,ThreeQuarterClade)
        
        ThreeQuartersimulatedtips<-simulatedtips_neutral[names(simulatedtips_neutral) %in% ThreeQuarterClade]
        
        ThreeQuartersimulatedtips<-droplevels(ThreeQuartersimulatedtips)
        
        selection_results[counter,1]<-tree_used
        selection_results[counter,2]<-4
        if(selection_variant<4) selection_results[counter,2]<-2
        selection_results[counter,3]<-selection_variant-3
        if(selection_variant<4) selection_results[counter,3]<-selection_variant
        selection_results[counter,4]<-repetition
        
        selection_results[counter,5]<-"ThreeQuarterSample"
        selection_results[counter,6]<-length(unique(ThreeQuartersimulatedtips))
        
        if (length(table(ThreeQuartersimulatedtips))!=1) { 
          
          
          numeric_simulatedtipds_neutral<-as.numeric(simulatedtips_neutral)
          names(numeric_simulatedtipds_neutral)<-names(simulatedtips_neutral)
          resultphylosiglambda<-phylosig(currenttree,numeric_simulatedtipds_neutral,method="lambda") 
          selection_results[counter,12]<-resultphylosiglambda$lambda
          resultphylosigK<-phylosig(currenttree,numeric_simulatedtipds_neutral,method="K") 
          selection_results[counter,13]<-resultphylosigK[1]
          
          #For the simulations where there are only two variants, we check whether by chance the variants
          #end up isolated in Clade A or associated with societies having different ecologies
          #For this we use Bayestraits Discrete, assessing the likelihood of independent versus dependent evolution
          if(length(table(ThreeQuartersimulatedtips))==2) {
            
            bayestraitsdiscretedata<-matrix(NA,nrow=length(ThreeQuartersimulatedtips),ncol=3)
            colnames(bayestraitsdiscretedata)<-c("Species","Simulated","CladeA")
            rownames(bayestraitsdiscretedata)<-c(1:length(ThreeQuartersimulatedtips))
            bayestraitsdiscretedata[,1]<-names(ThreeQuartersimulatedtips)
            bayestraitsdiscretedata[,2]<-as.numeric(as.character(ThreeQuartersimulatedtips))
            bayestraitsdiscretedata[,2]<-as.numeric(bayestraitsdiscretedata[,2])-as.numeric(min(bayestraitsdiscretedata[,2]))
            bayestraitsdiscretedata[bayestraitsdiscretedata[,2]!=0,2]<-1
            bayestraitsdiscretedata[,3]<- as.numeric(as.character(cladeAmember[names(cladeAmember) %in% names(ThreeQuartersimulatedtips)]))
            bayestraitsdiscretedata<-as.data.frame(bayestraitsdiscretedata)
            
            command_vec1 <- c("2", "1") #option 1 = 2 discrete independent; option 2 = 1 maximum likelihood
            results_1 <- bayestraits(bayestraitsdiscretedata, ThreeQuartertree, command_vec1)
            log_1 <- results_1$Log
            selection_results[counter,14]<-log_1$results$Lh
            
            command_vec2 <- c("3", "1") #option 1 = 3 discrete dependent; option 2 = 1 maximum likelihood
            results_2 <- bayestraits(bayestraitsdiscretedata, ThreeQuartertree, command_vec2)
            log_2 <- results_2$Log
            selection_results[counter,15]<-log_2$results$Lh
            
            
            bayestraitsdiscretedata<-matrix(NA,nrow=length(ThreeQuartersimulatedtips),ncol=3)
            colnames(bayestraitsdiscretedata)<-c("Species","Simulated","Ecologymember")
            rownames(bayestraitsdiscretedata)<-c(1:length(ThreeQuartersimulatedtips))
            bayestraitsdiscretedata[,1]<-names(ThreeQuartersimulatedtips)
            bayestraitsdiscretedata[,2]<-as.numeric(as.character(ThreeQuartersimulatedtips))
            bayestraitsdiscretedata[,2]<-as.numeric(bayestraitsdiscretedata[,2])-as.numeric(min(bayestraitsdiscretedata[,2]))
            bayestraitsdiscretedata[bayestraitsdiscretedata[,2]!=0,2]<-1
            bayestraitsdiscretedata[,3]<- as.numeric(as.character(Ecologymember[names(Ecologymember) %in% names(ThreeQuartersimulatedtips)]))
            bayestraitsdiscretedata[,3]<-as.numeric(as.character(bayestraitsdiscretedata[,3]))
            bayestraitsdiscretedata<-as.data.frame(bayestraitsdiscretedata)
            bayestraitsdiscretedata[is.na(bayestraitsdiscretedata[,3]),3]<-"-"
            command_vec3 <- c("2", "1") #option 1 = 2 discrete independent; option 2 = 1 maximum likelihood
            results_3 <- bayestraits(bayestraitsdiscretedata, ThreeQuartertree, command_vec3)
            log_3 <- results_3$Log
            selection_results[counter,16]<-log_3$results$Lh
            
            command_vec4 <- c("3", "1") #option 1 = 3 discrete dependent; option 2 = 1 maximum likelihood
            results_4 <- bayestraits(bayestraitsdiscretedata, ThreeQuartertree, command_vec4)
            log_4 <- results_4$Log
            selection_results[counter,17]<-log_4$results$Lh
            
          } #end of bayestraits 2 option
          
          
          
          
          
          
        }
        
        counter<-counter+1
        
        
        
        
        
        
        
        #-#-#-#-#-#-#-#-#-#-#-#
        
        # Subset the tree and the data to assume that 50% of data is randomly missing
        
        
        
        HalfClade <- currenttree$tip.label
        HalfClade <- sample(HalfClade)
        HalfClade <- HalfClade[1:round(nrow(data)/2,0)]
        
        Halftree<-keep.tip(currenttree,HalfClade)
        
        Halfsimulatedtips<-simulatedtips_neutral[names(simulatedtips_neutral) %in% HalfClade]
        
        Halfsimulatedtips<-droplevels(Halfsimulatedtips)
        
        selection_results[counter,1]<-tree_used
        selection_results[counter,2]<-4
        if(selection_variant<4) selection_results[counter,2]<-2
        selection_results[counter,3]<-selection_variant-3
        if(selection_variant<4) selection_results[counter,3]<-selection_variant          
        selection_results[counter,4]<-repetition
        
        selection_results[counter,5]<-"HalfSample"
        selection_results[counter,6]<-length(unique(Halfsimulatedtips))
        
        if (length(table(Halfsimulatedtips))!=1) { 
          
             
          numeric_simulatedtipds_neutral<-as.numeric(Halfsimulatedtips)
          names(numeric_simulatedtipds_neutral)<-names(Halfsimulatedtips)
          resultphylosiglambda<-phylosig(Halftree,numeric_simulatedtipds_neutral,method="lambda") 
          selection_results[counter,12]<-resultphylosiglambda$lambda
          resultphylosigK<-phylosig(Halftree,numeric_simulatedtipds_neutral,method="K") 
          selection_results[counter,13]<-resultphylosigK[1]
          
          
          #For the simulations where there are only two variants, we check whether by chance the variants
          #end up isolated in Clade A or associated with societies having different ecologies
          #For this we use Bayestraits Discrete, assessing the likelihood of independent versus dependent evolution
          if(length(table(Halfsimulatedtips))==2) {
            
            bayestraitsdiscretedata<-matrix(NA,nrow=length(Halfsimulatedtips),ncol=3)
            colnames(bayestraitsdiscretedata)<-c("Species","Simulated","CladeA")
            rownames(bayestraitsdiscretedata)<-c(1:length(Halfsimulatedtips))
            bayestraitsdiscretedata[,1]<-names(Halfsimulatedtips)
            bayestraitsdiscretedata[,2]<-as.numeric(as.character(Halfsimulatedtips))
            bayestraitsdiscretedata[,2]<-as.numeric(bayestraitsdiscretedata[,2])-as.numeric(min(bayestraitsdiscretedata[,2]))
            bayestraitsdiscretedata[bayestraitsdiscretedata[,2]!=0,2]<-1
            bayestraitsdiscretedata[,3]<- as.numeric(as.character(cladeAmember[names(cladeAmember) %in% names(Halfsimulatedtips)]))
            bayestraitsdiscretedata<-as.data.frame(bayestraitsdiscretedata)
            
            command_vec1 <- c("2", "1") #option 1 = 2 discrete independent; option 2 = 1 maximum likelihood
            results_1 <- bayestraits(bayestraitsdiscretedata, Halftree, command_vec1)
            log_1 <- results_1$Log
            selection_results[counter,14]<-log_1$results$Lh
            
            command_vec2 <- c("3", "1") #option 1 = 3 discrete dependent; option 2 = 1 maximum likelihood
            results_2 <- bayestraits(bayestraitsdiscretedata, Halftree, command_vec2)
            log_2 <- results_2$Log
            selection_results[counter,15]<-log_2$results$Lh
            
            
            bayestraitsdiscretedata<-matrix(NA,nrow=length(Halfsimulatedtips),ncol=3)
            colnames(bayestraitsdiscretedata)<-c("Species","Simulated","Ecologymember")
            rownames(bayestraitsdiscretedata)<-c(1:length(Halfsimulatedtips))
            bayestraitsdiscretedata[,1]<-names(Halfsimulatedtips)
            bayestraitsdiscretedata[,2]<-as.numeric(as.character(Halfsimulatedtips))
            bayestraitsdiscretedata[,2]<-as.numeric(bayestraitsdiscretedata[,2])-as.numeric(min(bayestraitsdiscretedata[,2]))
            bayestraitsdiscretedata[bayestraitsdiscretedata[,2]!=0,2]<-1
            bayestraitsdiscretedata[,3]<- as.numeric(as.character(Ecologymember[names(Ecologymember) %in% names(Halfsimulatedtips)]))
            bayestraitsdiscretedata[,3]<-as.numeric(as.character(bayestraitsdiscretedata[,3]))
            bayestraitsdiscretedata<-as.data.frame(bayestraitsdiscretedata)
            bayestraitsdiscretedata[is.na(bayestraitsdiscretedata[,3]),3]<-"-"
            command_vec3 <- c("2", "1") #option 1 = 2 discrete independent; option 2 = 1 maximum likelihood
            results_3 <- bayestraits(bayestraitsdiscretedata, Halftree, command_vec3)
            log_3 <- results_3$Log
            selection_results[counter,16]<-log_3$results$Lh
            
            command_vec4 <- c("3", "1") #option 1 = 3 discrete dependent; option 2 = 1 maximum likelihood
            results_4 <- bayestraits(bayestraitsdiscretedata, Halftree, command_vec4)
            log_4 <- results_4$Log
            selection_results[counter,17]<-log_4$results$Lh
            
          } #end of bayestraits 2 option
          
          
        }
        
        counter<-counter+1
        
        
        
        
        
        
        #-#-#-#-#-#-#-#-#-#-#-#
        
        # State dependent loss of samples - 50% from population of least frequent variant
        
        
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
        
        selection_results[counter,1]<-tree_used
        selection_results[counter,2]<-4
        if(selection_variant<4) selection_results[counter,2]<-2
        selection_results[counter,3]<-selection_variant-3
        if(selection_variant<4) selection_results[counter,3]<-selection_variant
        selection_results[counter,4]<-repetition
        
        selection_results[counter,5]<-"Rarelost"
        selection_results[counter,6]<-length(unique(Rarelostsimulatedtips))
        
        if (length(table(Rarelostsimulatedtips))!=1) { 
          
          # the phylogenetic reconstruction with the drift model, the model under which the data were simulated
            
          
          numeric_simulatedtipds_neutral<-as.numeric(Rarelostsimulatedtips)
          names(numeric_simulatedtipds_neutral)<-names(Rarelostsimulatedtips)
          resultphylosiglambda<-phylosig(Rarelosttree,numeric_simulatedtipds_neutral,method="lambda") 
          selection_results[counter,12]<-resultphylosiglambda$lambda
          resultphylosigK<-phylosig(Rarelosttree,numeric_simulatedtipds_neutral,method="K") 
          selection_results[counter,13]<-resultphylosigK[1]
          
          #For the simulations where there are only two variants, we check whether by chance the variants
          #end up isolated in Clade A or associated with societies having different ecologies
          #For this we use Bayestraits Discrete, assessing the likelihood of independent versus dependent evolution
          if(length(table(Rarelostsimulatedtips))==2) {
            
            bayestraitsdiscretedata<-matrix(NA,nrow=length(Rarelostsimulatedtips),ncol=3)
            colnames(bayestraitsdiscretedata)<-c("Species","Simulated","CladeA")
            rownames(bayestraitsdiscretedata)<-c(1:length(Rarelostsimulatedtips))
            bayestraitsdiscretedata[,1]<-names(Rarelostsimulatedtips)
            bayestraitsdiscretedata[,2]<-as.numeric(as.character(Rarelostsimulatedtips))
            bayestraitsdiscretedata[,2]<-as.numeric(bayestraitsdiscretedata[,2])-as.numeric(min(bayestraitsdiscretedata[,2]))
            bayestraitsdiscretedata[bayestraitsdiscretedata[,2]!=0,2]<-1
            bayestraitsdiscretedata[,3]<- as.numeric(as.character(cladeAmember[names(cladeAmember) %in% names(Rarelostsimulatedtips)]))
            bayestraitsdiscretedata<-as.data.frame(bayestraitsdiscretedata)
            
            command_vec1 <- c("2", "1") #option 1 = 2 discrete independent; option 2 = 1 maximum likelihood
            results_1 <- bayestraits(bayestraitsdiscretedata, Rarelosttree, command_vec1)
            log_1 <- results_1$Log
            selection_results[counter,14]<-log_1$results$Lh
            
            command_vec2 <- c("3", "1") #option 1 = 3 discrete dependent; option 2 = 1 maximum likelihood
            results_2 <- bayestraits(bayestraitsdiscretedata, Rarelosttree, command_vec2)
            log_2 <- results_2$Log
            selection_results[counter,15]<-log_2$results$Lh
            
            
            bayestraitsdiscretedata<-matrix(NA,nrow=length(Rarelostsimulatedtips),ncol=3)
            colnames(bayestraitsdiscretedata)<-c("Species","Simulated","Ecologymember")
            rownames(bayestraitsdiscretedata)<-c(1:length(Rarelostsimulatedtips))
            bayestraitsdiscretedata[,1]<-names(Rarelostsimulatedtips)
            bayestraitsdiscretedata[,2]<-as.numeric(as.character(Rarelostsimulatedtips))
            bayestraitsdiscretedata[,2]<-as.numeric(bayestraitsdiscretedata[,2])-as.numeric(min(bayestraitsdiscretedata[,2]))
            bayestraitsdiscretedata[bayestraitsdiscretedata[,2]!=0,2]<-1
            bayestraitsdiscretedata[,3]<- as.numeric(as.character(Ecologymember[names(Ecologymember) %in% names(Rarelostsimulatedtips)]))
            bayestraitsdiscretedata[,3]<-as.numeric(as.character(bayestraitsdiscretedata[,3]))
            bayestraitsdiscretedata<-as.data.frame(bayestraitsdiscretedata)
            bayestraitsdiscretedata[is.na(bayestraitsdiscretedata[,3]),3]<-"-"
            command_vec3 <- c("2", "1") #option 1 = 2 discrete independent; option 2 = 1 maximum likelihood
            results_3 <- bayestraits(bayestraitsdiscretedata, Rarelosttree, command_vec3)
            log_3 <- results_3$Log
            selection_results[counter,16]<-log_3$results$Lh
            
            command_vec4 <- c("3", "1") #option 1 = 3 discrete dependent; option 2 = 1 maximum likelihood
            results_4 <- bayestraits(bayestraitsdiscretedata, Rarelosttree, command_vec4)
            log_4 <- results_4$Log
            selection_results[counter,17]<-log_4$results$Lh
            
          } #end of bayestraits 2 option
          
          
          
          
        }
        
        counter<-counter+1
        
        
        
        
        
        
        #-#-#-#-#-#-#-#-#-#-#-#
        
        # State dependent loss of samples - 50% from population of most frequent variant
        
        
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
        
        selection_results[counter,1]<-tree_used
        selection_results[counter,2]<-4
        if(selection_variant<4) selection_results[counter,2]<-2
        selection_results[counter,3]<-selection_variant-3
        if(selection_variant<4) selection_results[counter,3]<-selection_variant
        selection_results[counter,4]<-repetition
        
        selection_results[counter,5]<-"Frequentlost"
        selection_results[counter,6]<-length(unique(Frequentlostsimulatedtips))
        
        if (length(table(Frequentlostsimulatedtips))!=1) { 
          
             
          
          numeric_simulatedtipds_neutral<-as.numeric(Frequentlostsimulatedtips)
          names(numeric_simulatedtipds_neutral)<-names(Frequentlostsimulatedtips)
          resultphylosiglambda<-phylosig(Frequentlosttree,numeric_simulatedtipds_neutral,method="lambda") 
          selection_results[counter,12]<-resultphylosiglambda$lambda
          resultphylosigK<-phylosig(Frequentlosttree,numeric_simulatedtipds_neutral,method="K") 
          selection_results[counter,13]<-resultphylosigK[1]
          
          #For the simulations where there are only two variants, we check whether by chance the variants
          #end up isolated in Clade A or associated with societies having different ecologies
          #For this we use Bayestraits Discrete, assessing the likelihood of independent versus dependent evolution
          if(length(table(Frequentlostsimulatedtips))==2) {
            
            bayestraitsdiscretedata<-matrix(NA,nrow=length(Frequentlostsimulatedtips),ncol=3)
            colnames(bayestraitsdiscretedata)<-c("Species","Simulated","CladeA")
            rownames(bayestraitsdiscretedata)<-c(1:length(Frequentlostsimulatedtips))
            bayestraitsdiscretedata[,1]<-names(Frequentlostsimulatedtips)
            bayestraitsdiscretedata[,2]<-as.numeric(as.character(Frequentlostsimulatedtips))
            bayestraitsdiscretedata[,2]<-as.numeric(bayestraitsdiscretedata[,2])-as.numeric(min(bayestraitsdiscretedata[,2]))
            bayestraitsdiscretedata[bayestraitsdiscretedata[,2]!=0,2]<-1
            bayestraitsdiscretedata[,3]<- as.numeric(as.character(cladeAmember[names(cladeAmember) %in% names(Frequentlostsimulatedtips)]))
            bayestraitsdiscretedata<-as.data.frame(bayestraitsdiscretedata)
            
            command_vec1 <- c("2", "1") #option 1 = 2 discrete independent; option 2 = 1 maximum likelihood
            results_1 <- bayestraits(bayestraitsdiscretedata, Frequentlosttree, command_vec1)
            log_1 <- results_1$Log
            selection_results[counter,14]<-log_1$results$Lh
            
            command_vec2 <- c("3", "1") #option 1 = 3 discrete dependent; option 2 = 1 maximum likelihood
            results_2 <- bayestraits(bayestraitsdiscretedata, Frequentlosttree, command_vec2)
            log_2 <- results_2$Log
            selection_results[counter,15]<-log_2$results$Lh
            
            
            bayestraitsdiscretedata<-matrix(NA,nrow=length(Frequentlostsimulatedtips),ncol=3)
            colnames(bayestraitsdiscretedata)<-c("Species","Simulated","Ecologymember")
            rownames(bayestraitsdiscretedata)<-c(1:length(Frequentlostsimulatedtips))
            bayestraitsdiscretedata[,1]<-names(Frequentlostsimulatedtips)
            bayestraitsdiscretedata[,2]<-as.numeric(as.character(Frequentlostsimulatedtips))
            bayestraitsdiscretedata[,2]<-as.numeric(bayestraitsdiscretedata[,2])-as.numeric(min(bayestraitsdiscretedata[,2]))
            bayestraitsdiscretedata[bayestraitsdiscretedata[,2]!=0,2]<-1
            bayestraitsdiscretedata[,3]<- as.numeric(as.character(Ecologymember[names(Ecologymember) %in% names(Frequentlostsimulatedtips)]))
            bayestraitsdiscretedata[,3]<-as.numeric(as.character(bayestraitsdiscretedata[,3]))
            bayestraitsdiscretedata<-as.data.frame(bayestraitsdiscretedata)
            bayestraitsdiscretedata[is.na(bayestraitsdiscretedata[,3]),3]<-"-"
            command_vec3 <- c("2", "1") #option 1 = 2 discrete independent; option 2 = 1 maximum likelihood
            results_3 <- bayestraits(bayestraitsdiscretedata, Frequentlosttree, command_vec3)
            log_3 <- results_3$Log
            selection_results[counter,16]<-log_3$results$Lh
            
            command_vec4 <- c("3", "1") #option 1 = 3 discrete dependent; option 2 = 1 maximum likelihood
            results_4 <- bayestraits(bayestraitsdiscretedata, Frequentlosttree, command_vec4)
            log_4 <- results_4$Log
            selection_results[counter,17]<-log_4$results$Lh
            
          } #end of bayestraits 2 option
          
          
          
          
        }
        
        counter<-counter+1
        
        
        
        
        
        
        
        # Loss according to ecoregion (50% of those living in forests or hunting/gathering)
        
        
        # Subset the tree and the data to assume that 50% of data is randomly missing
        
        
        
        
        randomdeletion<-function(x) (x*rbinom(1,1,0.5))
        Ecologytips<-sapply(as.numeric(Ecologymember),randomdeletion)
        
        
        if(Option=="WNAI") Ecologylosstips<-simulatedtips_neutral[Ecologytips==0]
        if(Option=="PamaNyungan") Ecologylosstips<-simulatedtips_neutral[Ecologytips==1]
        
        Ecologylosttree<-keep.tip(currenttree,names(Ecologylosstips))
        
        Ecologylosstips<-droplevels(Ecologylosstips)
        
        selection_results[counter,1]<-tree_used
        selection_results[counter,2]<-4
        if(selection_variant<4) selection_results[counter,2]<-2
        selection_results[counter,3]<-selection_variant-3
        if(selection_variant<4) selection_results[counter,3]<-selection_variant
        selection_results[counter,4]<-repetition
        
        selection_results[counter,5]<-"Ecologyloss"
        selection_results[counter,6]<-length(unique(Ecologylosstips))
        
        if (length(table(Ecologylosstips))!=1) { 
          
         
          numeric_simulatedtipds_neutral<-as.numeric(Ecologylosstips)
          names(numeric_simulatedtipds_neutral)<-names(Ecologylosstips)
          resultphylosiglambda<-phylosig(Ecologylosttree,numeric_simulatedtipds_neutral,method="lambda") 
          selection_results[counter,12]<-resultphylosiglambda$lambda
          resultphylosigK<-phylosig(Ecologylosttree,numeric_simulatedtipds_neutral,method="K") 
          selection_results[counter,13]<-resultphylosigK[1]
          
          
          #For the simulations where there are only two variants, we check whether by chance the variants
          #end up isolated in Clade A or associated with societies having different ecologies
          #For this we use Bayestraits Discrete, assessing the likelihood of independent versus dependent evolution
          if(length(table(Ecologylosstips))==2) {
            
            bayestraitsdiscretedata<-matrix(NA,nrow=length(Ecologylosstips),ncol=3)
            colnames(bayestraitsdiscretedata)<-c("Species","Simulated","CladeA")
            rownames(bayestraitsdiscretedata)<-c(1:length(Ecologylosstips))
            bayestraitsdiscretedata[,1]<-names(Ecologylosstips)
            bayestraitsdiscretedata[,2]<-as.numeric(as.character(Ecologylosstips))
            bayestraitsdiscretedata[,2]<-as.numeric(bayestraitsdiscretedata[,2])-as.numeric(min(bayestraitsdiscretedata[,2]))
            bayestraitsdiscretedata[bayestraitsdiscretedata[,2]!=0,2]<-1
            bayestraitsdiscretedata[,3]<- as.numeric(as.character(cladeAmember[names(cladeAmember) %in% names(Ecologylosstips)]))
            bayestraitsdiscretedata<-as.data.frame(bayestraitsdiscretedata)
            
            command_vec1 <- c("2", "1") #option 1 = 2 discrete independent; option 2 = 1 maximum likelihood
            results_1 <- bayestraits(bayestraitsdiscretedata, Ecologylosttree, command_vec1)
            log_1 <- results_1$Log
            selection_results[counter,14]<-log_1$results$Lh
            
            command_vec2 <- c("3", "1") #option 1 = 3 discrete dependent; option 2 = 1 maximum likelihood
            results_2 <- bayestraits(bayestraitsdiscretedata, Ecologylosttree, command_vec2)
            log_2 <- results_2$Log
            selection_results[counter,15]<-log_2$results$Lh
            
            
            bayestraitsdiscretedata<-matrix(NA,nrow=length(Ecologylosstips),ncol=3)
            colnames(bayestraitsdiscretedata)<-c("Species","Simulated","Ecologymember")
            rownames(bayestraitsdiscretedata)<-c(1:length(Ecologylosstips))
            bayestraitsdiscretedata[,1]<-names(Ecologylosstips)
            bayestraitsdiscretedata[,2]<-as.numeric(as.character(Ecologylosstips))
            bayestraitsdiscretedata[,2]<-as.numeric(bayestraitsdiscretedata[,2])-as.numeric(min(bayestraitsdiscretedata[,2]))
            bayestraitsdiscretedata[bayestraitsdiscretedata[,2]!=0,2]<-1
            bayestraitsdiscretedata[,3]<- as.numeric(as.character(Ecologymember[names(Ecologymember) %in% names(Ecologylosstips)]))
            bayestraitsdiscretedata[,3]<-as.numeric(as.character(bayestraitsdiscretedata[,3]))
            bayestraitsdiscretedata<-as.data.frame(bayestraitsdiscretedata)
            bayestraitsdiscretedata[is.na(bayestraitsdiscretedata[,3]),3]<-"-"
            command_vec3 <- c("2", "1") #option 1 = 2 discrete independent; option 2 = 1 maximum likelihood
            results_3 <- bayestraits(bayestraitsdiscretedata, Ecologylosttree, command_vec3)
            log_3 <- results_3$Log
            selection_results[counter,16]<-log_3$results$Lh
            
            command_vec4 <- c("3", "1") #option 1 = 3 discrete dependent; option 2 = 1 maximum likelihood
            results_4 <- bayestraits(bayestraitsdiscretedata, Ecologylosttree, command_vec4)
            log_4 <- results_4$Log
            selection_results[counter,17]<-log_4$results$Lh
            
          } #end of bayestraits 2 option
          
          
          
        }
        
        counter<-counter+1
        
        
        
        
        
        
        
        
        #--#--#--#--#--#--#        #-#-#-#-#-#-#-#-#-#-#-#-#
        
        
        # Horizontal transfer, 10% of societies change, randomly adopt another state
        
        
        
        societiestochange<-sample(simulatedtips_neutral,round(0.1*nrow(data),0))
        societiestochange[societiestochange==1]<-NA
        societiestochange[societiestochange==2]<-1
        societiestochange[is.na(societiestochange)==T]<-2
        horizontal10tips<-simulatedtips_neutral
        horizontal10tips[names(horizontal10tips) %in% names(societiestochange)]<-societiestochange
        
        horizontal10tips<-droplevels(horizontal10tips)
        
        
        #Start filling in the data in the respective column
        selection_results[counter,1]<-tree_used
        selection_results[counter,2]<-4
        if(selection_variant<4) selection_results[counter,2]<-2
        selection_results[counter,3]<-selection_variant-3
        if(selection_variant<4) selection_results[counter,3]<-selection_variant
        
        selection_results[counter,4]<-repetition
        
        
        #Start with the full sample of  societies
        
        selection_results[counter,5]<-"horizontal10"
        selection_results[counter,6]<-length(unique(horizontal10tips))
        
        if (length(table(horizontal10tips))!=1) { 
          
              
          numeric_simulatedtipds_neutral<-as.numeric(horizontal10tips)
          names(numeric_simulatedtipds_neutral)<-names(horizontal10tips)
          resultphylosiglambda<-phylosig(currenttree,numeric_simulatedtipds_neutral,method="lambda") 
          selection_results[counter,12]<-resultphylosiglambda$lambda
          resultphylosigK<-phylosig(currenttree,numeric_simulatedtipds_neutral,method="K") 
          selection_results[counter,13]<-resultphylosigK[1]
          
          
          #For the simulations where there are only two variants, we check whether by chance the variants
          #end up isolated in Clade A or associated with societies living in different ecologies
          #For this we use Bayestraits Discrete, assessing the likelihood of independent versus dependent evolution
          if(length(table(horizontal10tips))==2) {
            
            bayestraitsdiscretedata<-matrix(NA,nrow=length(horizontal10tips),ncol=3)
            colnames(bayestraitsdiscretedata)<-c("Species","Simulated","CladeA")
            rownames(bayestraitsdiscretedata)<-c(1:length(horizontal10tips))
            bayestraitsdiscretedata[,1]<-names(horizontal10tips)
            bayestraitsdiscretedata[,2]<-as.numeric(as.character(horizontal10tips))
            bayestraitsdiscretedata[,2]<-as.numeric(bayestraitsdiscretedata[,2])-as.numeric(min(bayestraitsdiscretedata[,2]))
            bayestraitsdiscretedata[bayestraitsdiscretedata[,2]!=0,2]<-1
            bayestraitsdiscretedata[,3]<- as.numeric(as.character(cladeAmember[names(cladeAmember) %in% names(horizontal10tips)]))
            bayestraitsdiscretedata<-as.data.frame(bayestraitsdiscretedata)
            
            command_vec1 <- c("2", "1") #option 1 = 2 discrete independent; option 2 = 1 maximum likelihood
            results_1 <- bayestraits(bayestraitsdiscretedata, currenttree, command_vec1)
            log_1 <- results_1$Log
            selection_results[counter,14]<-log_1$results$Lh
            
            command_vec2 <- c("3", "1") #option 1 = 3 discrete dependent; option 2 = 1 maximum likelihood
            results_2 <- bayestraits(bayestraitsdiscretedata, currenttree, command_vec2)
            log_2 <- results_2$Log
            selection_results[counter,15]<-log_2$results$Lh
            
            
            bayestraitsdiscretedata<-matrix(NA,nrow=length(horizontal10tips),ncol=3)
            colnames(bayestraitsdiscretedata)<-c("Species","Simulated","Ecologymember")
            rownames(bayestraitsdiscretedata)<-c(1:length(horizontal10tips))
            bayestraitsdiscretedata[,1]<-names(horizontal10tips)
            bayestraitsdiscretedata[,2]<-as.numeric(as.character(horizontal10tips))
            bayestraitsdiscretedata[,2]<-as.numeric(bayestraitsdiscretedata[,2])-as.numeric(min(bayestraitsdiscretedata[,2]))
            bayestraitsdiscretedata[bayestraitsdiscretedata[,2]!=0,2]<-1
            bayestraitsdiscretedata[,3]<- as.numeric(as.character(Ecologymember[names(Ecologymember) %in% names(horizontal10tips)]))
            bayestraitsdiscretedata[,3]<-as.numeric(as.character(bayestraitsdiscretedata[,3]))
            bayestraitsdiscretedata<-as.data.frame(bayestraitsdiscretedata)
            bayestraitsdiscretedata[is.na(bayestraitsdiscretedata[,3]),3]<-"-"
            
            command_vec3 <- c("2", "1") #option 1 = 2 discrete independent; option 2 = 1 maximum likelihood
            results_3 <- bayestraits(bayestraitsdiscretedata, currenttree, command_vec3)
            log_3 <- results_3$Log
            selection_results[counter,16]<-log_3$results$Lh
            
            command_vec4 <- c("3", "1") #option 1 = 3 discrete dependent; option 2 = 1 maximum likelihood
            results_4 <- bayestraits(bayestraitsdiscretedata, currenttree, command_vec4)
            log_4 <- results_4$Log
            selection_results[counter,17]<-log_4$results$Lh
            
          } #end of bayestraits option
          
               
          
        }
        
        counter<-counter+1
        
        
        
        
        
        
        
        # Horizontal transfer, 10% of societies change, adopt state most frequent within their own clade
        
        cladeAsimulatedtips<-simulatedtips_neutral[names(simulatedtips_neutral) %in% cladeA]
        cladeAtree<-keep.tip(currenttree,cladeA)
        
        cladeAsimulatedtips<-droplevels(cladeAsimulatedtips)
        
        
        if (length(table(cladeAsimulatedtips))!=1) {
          
          cladeAtochange<-sample(cladeAsimulatedtips[cladeAsimulatedtips==names(table(cladeAsimulatedtips)[2])],min(round(0.1*nrow(data),0),table(cladeAsimulatedtips)[2]))
          levels(cladeAtochange)<-c("1","2","3","4")
          cladeAtochange[]<-1
          horizontalClade10tips<-simulatedtips_neutral
          horizontalClade10tips[names(horizontalClade10tips) %in% names(cladeAtochange)]<-cladeAtochange
          
          horizontalClade10tips<-droplevels(horizontalClade10tips)
          
          #Start filling in the data in the respective column
          selection_results[counter,1]<-tree_used
          selection_results[counter,2]<-4
          if(selection_variant<4) selection_results[counter,2]<-2
          selection_results[counter,3]<-selection_variant-3
          if(selection_variant<4) selection_results[counter,3]<-selection_variant
          selection_results[counter,4]<-repetition
          
          
          #Start with the full sample of societies
          
          selection_results[counter,5]<-"horizontalClade10"
          selection_results[counter,6]<-length(unique(horizontalClade10tips))
          
          if (length(table(horizontalClade10tips))!=1) { 
            
            
             
            numeric_simulatedtipds_neutral<-as.numeric(horizontalClade10tips)
            names(numeric_simulatedtipds_neutral)<-names(horizontalClade10tips)
            resultphylosiglambda<-phylosig(currenttree,numeric_simulatedtipds_neutral,method="lambda") 
            selection_results[counter,12]<-resultphylosiglambda$lambda
            resultphylosigK<-phylosig(currenttree,numeric_simulatedtipds_neutral,method="K") 
            selection_results[counter,13]<-resultphylosigK[1]
            
            
            #For the simulations where there are only two variants, we check whether by chance the variants
            #end up isolated in Clade A or associated with societies living in ecologies
            #For this we use Bayestraits Discrete, assessing the likelihood of independent versus dependent evolution
            if(length(table(horizontalClade10tips))==2) {
              
              bayestraitsdiscretedata<-matrix(NA,nrow=length(horizontalClade10tips),ncol=3)
              colnames(bayestraitsdiscretedata)<-c("Species","Simulated","CladeA")
              rownames(bayestraitsdiscretedata)<-c(1:length(horizontalClade10tips))
              bayestraitsdiscretedata[,1]<-names(horizontalClade10tips)
              bayestraitsdiscretedata[,2]<-as.numeric(as.character(horizontalClade10tips))
              bayestraitsdiscretedata[,2]<-as.numeric(bayestraitsdiscretedata[,2])-as.numeric(min(bayestraitsdiscretedata[,2]))
              bayestraitsdiscretedata[bayestraitsdiscretedata[,2]!=0,2]<-1
              bayestraitsdiscretedata[,3]<- as.numeric(as.character(cladeAmember[names(cladeAmember) %in% names(horizontalClade10tips)]))
              bayestraitsdiscretedata<-as.data.frame(bayestraitsdiscretedata)
              
              command_vec1 <- c("2", "1") #option 1 = 2 discrete independent; option 2 = 1 maximum likelihood
              results_1 <- bayestraits(bayestraitsdiscretedata, currenttree, command_vec1)
              log_1 <- results_1$Log
              selection_results[counter,14]<-log_1$results$Lh
              
              command_vec2 <- c("3", "1") #option 1 = 3 discrete dependent; option 2 = 1 maximum likelihood
              results_2 <- bayestraits(bayestraitsdiscretedata, currenttree, command_vec2)
              log_2 <- results_2$Log
              selection_results[counter,15]<-log_2$results$Lh
              
              
              bayestraitsdiscretedata<-matrix(NA,nrow=length(horizontalClade10tips),ncol=3)
              colnames(bayestraitsdiscretedata)<-c("Species","Simulated","Ecologymember")
              rownames(bayestraitsdiscretedata)<-c(1:length(horizontalClade10tips))
              bayestraitsdiscretedata[,1]<-names(horizontalClade10tips)
              bayestraitsdiscretedata[,2]<-as.numeric(as.character(horizontalClade10tips))
              bayestraitsdiscretedata[,2]<-as.numeric(bayestraitsdiscretedata[,2])-as.numeric(min(bayestraitsdiscretedata[,2]))
              bayestraitsdiscretedata[bayestraitsdiscretedata[,2]!=0,2]<-1
              bayestraitsdiscretedata[,3]<- as.numeric(as.character(Ecologymember[names(Ecologymember) %in% names(horizontalClade10tips)]))
              bayestraitsdiscretedata[,3]<-as.numeric(as.character(bayestraitsdiscretedata[,3]))
              bayestraitsdiscretedata<-as.data.frame(bayestraitsdiscretedata)
              bayestraitsdiscretedata[is.na(bayestraitsdiscretedata[,3]),3]<-"-"
              
              command_vec3 <- c("2", "1") #option 1 = 2 discrete independent; option 2 = 1 maximum likelihood
              results_3 <- bayestraits(bayestraitsdiscretedata, currenttree, command_vec3)
              log_3 <- results_3$Log
              selection_results[counter,16]<-log_3$results$Lh
              
              command_vec4 <- c("3", "1") #option 1 = 3 discrete dependent; option 2 = 1 maximum likelihood
              results_4 <- bayestraits(bayestraitsdiscretedata, currenttree, command_vec4)
              log_4 <- results_4$Log
              selection_results[counter,17]<-log_4$results$Lh
              
            } #end of bayestraits option
            
            
          
            
          } #end of checking whether there are more than 1 variants in the sample
          
        } #end of checking whether there are more than 1 variants in Clade A
        
        counter<-counter+1
        
        
        
        
        # Horizontal transfer, 10% of societies change, adopt state most frequent within their own ecoregion
        
        
        
        Ecologysimulatedtips<-simulatedtips_neutral[names(simulatedtips_neutral) %in% names(Ecologymember[Ecologymember==1])]
        
        Ecologytochange<-sample(Ecologysimulatedtips[Ecologysimulatedtips==names(table(Ecologysimulatedtips)[2])],min(17,table(Ecologysimulatedtips)[2]))
        levels(Ecologytochange)<-c("1","2","3","4")
        Ecologytochange[]<-1
        horizontalEcology10tips<-simulatedtips_neutral
        horizontalEcology10tips[names(horizontalEcology10tips) %in% names(Ecologytochange)]<-Ecologytochange
        
        horizontalEcology10tips<-droplevels(horizontalEcology10tips)
        
        #Start filling in the data in the respective column
        selection_results[counter,1]<-tree_used
        selection_results[counter,2]<-4
        if(selection_variant<4) selection_results[counter,2]<-2
        selection_results[counter,3]<-selection_variant-3
        if(selection_variant<4) selection_results[counter,3]<-selection_variant
        selection_results[counter,4]<-repetition
        
        
        #Start with the full sample 
        
        selection_results[counter,5]<-"horizontalEcology10"
        selection_results[counter,6]<-length(unique(horizontalEcology10tips))
        
        if (length(table(horizontalEcology10tips))!=1) { 
          
          
           
          numeric_simulatedtipds_neutral<-as.numeric(horizontalEcology10tips)
          names(numeric_simulatedtipds_neutral)<-names(horizontalEcology10tips)
          resultphylosiglambda<-phylosig(currenttree,numeric_simulatedtipds_neutral,method="lambda") 
          selection_results[counter,12]<-resultphylosiglambda$lambda
          resultphylosigK<-phylosig(currenttree,numeric_simulatedtipds_neutral,method="K") 
          selection_results[counter,13]<-resultphylosigK[1]
          
          
          #For the simulations where there are only two variants, we check whether by chance the variants
          #end up isolated in Clade A or associated with societies living in forest ecoregions
          #For this we use Bayestraits Discrete, assessing the likelihood of independent versus dependent evolution
          if(length(table(horizontalEcology10tips))==2) {
            
            bayestraitsdiscretedata<-matrix(NA,nrow=length(horizontalEcology10tips),ncol=3)
            colnames(bayestraitsdiscretedata)<-c("Species","Simulated","CladeA")
            rownames(bayestraitsdiscretedata)<-c(1:length(horizontalEcology10tips))
            bayestraitsdiscretedata[,1]<-names(horizontalEcology10tips)
            bayestraitsdiscretedata[,2]<-as.numeric(as.character(horizontalEcology10tips))
            bayestraitsdiscretedata[,2]<-as.numeric(bayestraitsdiscretedata[,2])-as.numeric(min(bayestraitsdiscretedata[,2]))
            bayestraitsdiscretedata[bayestraitsdiscretedata[,2]!=0,2]<-1
            bayestraitsdiscretedata[,3]<- as.numeric(as.character(cladeAmember[names(cladeAmember) %in% names(horizontalEcology10tips)]))
            bayestraitsdiscretedata<-as.data.frame(bayestraitsdiscretedata)
            
            command_vec1 <- c("2", "1") #option 1 = 2 discrete independent; option 2 = 1 maximum likelihood
            results_1 <- bayestraits(bayestraitsdiscretedata, currenttree, command_vec1)
            log_1 <- results_1$Log
            selection_results[counter,14]<-log_1$results$Lh
            
            command_vec2 <- c("3", "1") #option 1 = 3 discrete dependent; option 2 = 1 maximum likelihood
            results_2 <- bayestraits(bayestraitsdiscretedata, currenttree, command_vec2)
            log_2 <- results_2$Log
            selection_results[counter,15]<-log_2$results$Lh
            
            
            bayestraitsdiscretedata<-matrix(NA,nrow=length(horizontalEcology10tips),ncol=3)
            colnames(bayestraitsdiscretedata)<-c("Species","Simulated","Ecologymember")
            rownames(bayestraitsdiscretedata)<-c(1:length(horizontalEcology10tips))
            bayestraitsdiscretedata[,1]<-names(horizontalEcology10tips)
            bayestraitsdiscretedata[,2]<-as.numeric(as.character(horizontalEcology10tips))
            bayestraitsdiscretedata[,2]<-as.numeric(bayestraitsdiscretedata[,2])-as.numeric(min(bayestraitsdiscretedata[,2]))
            bayestraitsdiscretedata[bayestraitsdiscretedata[,2]!=0,2]<-1
            bayestraitsdiscretedata[,3]<- as.numeric(as.character(Ecologymember[names(Ecologymember) %in% names(horizontalEcology10tips)]))
            bayestraitsdiscretedata[,3]<-as.numeric(as.character(bayestraitsdiscretedata[,3]))
            bayestraitsdiscretedata<-as.data.frame(bayestraitsdiscretedata)
            bayestraitsdiscretedata[is.na(bayestraitsdiscretedata[,3]),3]<-"-"
            
            command_vec3 <- c("2", "1") #option 1 = 2 discrete independent; option 2 = 1 maximum likelihood
            results_3 <- bayestraits(bayestraitsdiscretedata, currenttree, command_vec3)
            log_3 <- results_3$Log
            selection_results[counter,16]<-log_3$results$Lh
            
            command_vec4 <- c("3", "1") #option 1 = 3 discrete dependent; option 2 = 1 maximum likelihood
            results_4 <- bayestraits(bayestraitsdiscretedata, currenttree, command_vec4)
            log_4 <- results_4$Log
            selection_results[counter,17]<-log_4$results$Lh
            
          } #end of bayestraits option
          
          
         
          
        }
        
        counter<-counter+1
        
        
        
        
        
        # Horizontal transfer, 10% of societies change, adopt state of a close neighbor
        
        
        
        Latitudedistance<-as.matrix(dist(locationdata$Latitude))
        Longitudedistance<-as.matrix(dist(locationdata$Longitude))
        Totaldistance<-Latitudedistance+Longitudedistance
        diag(Totaldistance)<-1000
        colnames(Totaldistance)<-row.names(locationdata)
        row.names(Totaldistance)<-row.names(locationdata)
        subsettochange<-sample(simulatedtips_neutral,round(0.15*nrow(data),0))
        
        for (findtheneighbor in 1:round(0.15*nrow(data),0)) {
          
          findtherow<-Totaldistance[,colnames(Totaldistance)==names(subsettochange[findtheneighbor])]==min(Totaldistance[,colnames(Totaldistance)==names(subsettochange[findtheneighbor]) ] ,na.rm=T)
          
          
          subsettochange[findtheneighbor]<-simulatedtips_neutral[names(simulatedtips_neutral)==min(names(findtherow[findtherow==TRUE]))]
          
        }
        
        horizontalNeighbor10tips<-simulatedtips_neutral
        horizontalNeighbor10tips[names(horizontalNeighbor10tips) %in% names(subsettochange)]<-subsettochange
        
        horizontalNeighbor10tips<-droplevels(horizontalNeighbor10tips)
        
        #Start filling in the data in the respective column
        selection_results[counter,1]<-tree_used
        selection_results[counter,2]<-4
        if(selection_variant<4) selection_results[counter,2]<-2
        selection_results[counter,3]<-selection_variant-3
        if(selection_variant<4) selection_results[counter,3]<-selection_variant
        selection_results[counter,4]<-repetition
        
        
        #Start with the full sample of societies
        
        selection_results[counter,5]<-"horizontalNeighbor10"
        selection_results[counter,6]<-length(unique(horizontalNeighbor10tips))
        
        if (length(table(horizontalNeighbor10tips))!=1) { 
          
          numeric_simulatedtipds_neutral<-as.numeric(horizontalNeighbor10tips)
          names(numeric_simulatedtipds_neutral)<-names(horizontalNeighbor10tips)
          resultphylosiglambda<-phylosig(currenttree,numeric_simulatedtipds_neutral,method="lambda") 
          selection_results[counter,12]<-resultphylosiglambda$lambda
          resultphylosigK<-phylosig(currenttree,numeric_simulatedtipds_neutral,method="K") 
          selection_results[counter,13]<-resultphylosigK[1]
          
          
          #For the simulations where there are only two variants, we check whether by chance the variants
          #end up isolated in Clade A or associated with societies living in different ecologies
          #For this we use Bayestraits Discrete, assessing the likelihood of independent versus dependent evolution
          if(length(table(horizontalNeighbor10tips))==2) {
            
            bayestraitsdiscretedata<-matrix(NA,nrow=length(horizontalNeighbor10tips),ncol=3)
            colnames(bayestraitsdiscretedata)<-c("Species","Simulated","CladeA")
            rownames(bayestraitsdiscretedata)<-c(1:length(horizontalNeighbor10tips))
            bayestraitsdiscretedata[,1]<-names(horizontalNeighbor10tips)
            bayestraitsdiscretedata[,2]<-as.numeric(as.character(horizontalNeighbor10tips))
            bayestraitsdiscretedata[,2]<-as.numeric(bayestraitsdiscretedata[,2])-as.numeric(min(bayestraitsdiscretedata[,2]))
            bayestraitsdiscretedata[bayestraitsdiscretedata[,2]!=0,2]<-1
            bayestraitsdiscretedata[,3]<- as.numeric(as.character(cladeAmember[names(cladeAmember) %in% names(horizontalNeighbor10tips)]))
            bayestraitsdiscretedata<-as.data.frame(bayestraitsdiscretedata)
            
            command_vec1 <- c("2", "1") #option 1 = 2 discrete independent; option 2 = 1 maximum likelihood
            results_1 <- bayestraits(bayestraitsdiscretedata, currenttree, command_vec1)
            log_1 <- results_1$Log
            selection_results[counter,14]<-log_1$results$Lh
            
            command_vec2 <- c("3", "1") #option 1 = 3 discrete dependent; option 2 = 1 maximum likelihood
            results_2 <- bayestraits(bayestraitsdiscretedata, currenttree, command_vec2)
            log_2 <- results_2$Log
            selection_results[counter,15]<-log_2$results$Lh
            
            
            bayestraitsdiscretedata<-matrix(NA,nrow=length(horizontalNeighbor10tips),ncol=3)
            colnames(bayestraitsdiscretedata)<-c("Species","Simulated","Ecologymember")
            rownames(bayestraitsdiscretedata)<-c(1:length(horizontalNeighbor10tips))
            bayestraitsdiscretedata[,1]<-names(horizontalNeighbor10tips)
            bayestraitsdiscretedata[,2]<-as.numeric(as.character(horizontalNeighbor10tips))
            bayestraitsdiscretedata[,2]<-as.numeric(bayestraitsdiscretedata[,2])-as.numeric(min(bayestraitsdiscretedata[,2]))
            bayestraitsdiscretedata[bayestraitsdiscretedata[,2]!=0,2]<-1
            bayestraitsdiscretedata[,3]<- as.numeric(as.character(Ecologymember[names(Ecologymember) %in% names(horizontalNeighbor10tips)]))
            bayestraitsdiscretedata[,3]<-as.numeric(as.character(bayestraitsdiscretedata[,3]))
            bayestraitsdiscretedata<-as.data.frame(bayestraitsdiscretedata)
            bayestraitsdiscretedata[is.na(bayestraitsdiscretedata[,3]),3]<-"-"
            
            command_vec3 <- c("2", "1") #option 1 = 2 discrete independent; option 2 = 1 maximum likelihood
            results_3 <- bayestraits(bayestraitsdiscretedata, currenttree, command_vec3)
            log_3 <- results_3$Log
            selection_results[counter,16]<-log_3$results$Lh
            
            command_vec4 <- c("3", "1") #option 1 = 3 discrete dependent; option 2 = 1 maximum likelihood
            results_4 <- bayestraits(bayestraitsdiscretedata, currenttree, command_vec4)
            log_4 <- results_4$Log
            selection_results[counter,17]<-log_4$results$Lh
            
          } #end of bayestraits option
          
          
       
          
        }
        
        counter<-counter+1
        
        
        print("runningcounter")
        print(counter)
        print("repetition")
        print(repetition)
        print("currenttree")
        print(tree_used)
        print("selection_variant")
        print(selection_variant)
        
        
      } #only need to run if there is more than one variant  
      
    } # end of the repetition loop
    
    if(Option=="WNAI") {
      if(condition_variant==1) {write.csv(selection_results[1:counter,],file="multipleselectionresults_WNAI_clade.csv") }
      if(condition_variant==2) {write.csv(selection_results[1:counter,],file="multipleselectionresults_WNAI_ecology.csv") }
      }
    
    if(Option=="PamaNyungan") {
      if(condition_variant==1) {write.csv(selection_results[1:counter,],file="multipleselectionresults_PamaNyungan_clade.csv") }
      if(condition_variant==2) {write.csv(selection_results[1:counter,],file="multipleselectionresults_PamaNyungan_ecology.csv") }
      }
    
  }  # end of the drift variants model
  
  
} #end of the tree variant loop
  
} #end of the condition variant loop

options(warn = oldw)

#In case you want to save the results again
if(Option=="WNAI") {write.csv(selection_results[1:counter,],file="SimulateCulture_MultipleSelectionResults_WNAI_clade.csv") }
if(Option=="PamaNyungan") {write.csv(selection_results[1:counter,],file="SimulateCulture_MultipleSelectionResults_PamaNyungan_clade.csv") }

