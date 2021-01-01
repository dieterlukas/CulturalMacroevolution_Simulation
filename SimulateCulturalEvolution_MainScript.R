

#------------------------------------------------------------------------------------------
# Simulations of the evolution of cultural traits across phylogenies
#------------------------------------------------------------------------------------------

# These simulations and their outcomes are described in:
# Lukas, D., Towner, M., & Mulder, M. B. (2020). The Potential to Infer the Historical Pattern of Cultural Macroevolution as Illustrated by the Western North American Indian Societies.
# https://osf.io/preprints/socarxiv/tjvgy/



#------------------------------------------------------------------------------------------
# This set of simulations and analyses aims to assess the potential error rates
# of phylogenetic reconstructions of the historical evolution of cultural traits

# The code simulates changes in an arbitrary cultural trait across a known phylogeny. The trait can be
# either continuous or discrete, changing randomly or differently across the phylogeny.
# The resulting distribution of the different states of the trait are analysed with a variety of methods
# that are frequently used in biological and cultural phylogenetics. A false positive error occurs
# if such an analysis wrongly infers that the trait has been changing in one particular direction or
# that changes in the trait occurred at different rates in different parts of the tree even though the simulation
# did not specify any direction or differences across the tree. A false negative error occurs if such an 
# analysis wrongly infers that the trait has been changing equally in all directions throughout the tree 
# even though the simulation specified differences and changes in particular directions.

# The analyses examine the potential false errors linked to three inferences:
# 1) wrongly identifying or missing directional changes of the trait:
#         Type 1 error / false positive: an analysis wrongly supports the inference that changes to or from one state or in
#         one direction have occurred significantly more likely than other changes - even though the simulation is set up
#         such that changes between all states are equally likely.
#         Type 2 error / false negative: an analysis wrongly supports the inferences that all changes have occurred
#         equally likely - even though the simulation is set up such that some changes occur more frequently than others.
# 2) wrongly identifying or missing lineage differences in changes in the trait:
#         Type 1 error / false positive: an analysis wrongly supports the inference that some changes
#         occurred more or less frequently in one part of the tree compared to the rest of the tree - even
#         though the simulation is set up such that changes occur at the same rate throughout the tree
#         Type 2 error / false negative: an analysis wrongly supports the inferences that changes did not differ
#         between different parts of the tree - even though the simulation is set up such that some changes are
#         more frequently in one clade than the rest of the tree
# 3) wronlgy identifying or missing ecological associations of the distribution of the trait:
#         Type 1 error / false positive: an analysis wrongly supports that changes to or from one state have occurred
#         significantly more often in lineages that have a particular ecology - even though the simulation
#         is set up such that changes in the trait are not associated with ecological differences
#         Type 2 error / false negative: an analysis wrongly supports that changes occurred equally throughout
#         the tree - even though the simulation is set up such that some changes occur more in lineages
#         associated with specific ecological conditions.


# The setup includes further modifications to assess whether the error rates are influenced by
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


#------------------------------------------------------------------------------------------


# First, please select a folder on your computer as the working directory (setwd or Session -> Set Working Directory)
# All output will be saved to this folder. In addition, you will need to have a copy of the software Bayestraits
# in this folder for the simulations of the discrete traits.

# The simulation relies on a number of different packages in R.

#Load necessary libraries
# The code relies on a number of phylogenetic packages to manipulate and analyse the data
library(ape)
library(geiger)
library(phytools)
library(OUwie)
library(caper)

# The code also uses some packages that facilitate loading and structuring the data
library(dplyr)
library(readr)

# Finally, the simulation relies on the software Bayestraits (http://www.evolution.rdg.ac.uk/SoftwareMain.html)
# It starts the software from within R using the package Bayestraits Wrapper
# To install the package, follow the instructions by the developer (http://www.randigriffin.com/projects/btw.html)
# It is also necessary that you have Bayestraits on your computer, and that you set the working directory in R to a location with a copy of Bayestraits
library(btw)



#------------------------------------------------------------------------------------------
# In our manuscript, we illustrate this simulation approache with two sets of societies here as example:
# Western North American Indigeneous societies (WNAI) and societies from Australia who speak a language of the Pama Nyungan family

# To run the simulations across a phylogeny of Western North American Indigeneous societies, run this next line
Option<-"WNAI"

# To run the simulations across a phylogeny of Pama PamaNyungan speaking societies, run this next line
Option<-"PamaNyungan"



# The code can be adapted to assess the potential risk of false positive errors in a user-provided dataset.
# For this, you need a phylogeny in nexus format, and, if relevant, information on an ecological variable
# for the societies in the sample. Load these instead of the example files listed below.

# The code is provided in separate files in a folder on GitHub. They are being called directly from 
# this master file. The additional code files contains further comments and explanations. They can 
# be used directly to simulate cultural traits on user-provided datasets.

#------------------------------------------------------------------------------------------


# All data necessary for these examples are available on GitHub, and they can be directly loaded here.
# This includes the phylogeny, the location of each society, and the coding of an ecological predictor variable

source("https://raw.githubusercontent.com/dieterlukas/CulturalMacroevolution_Simulation/master/Code/SimulateCulturalEvolution_LoadFiles.R")

#------------------------------------------------------------------------------------------

# To determine whether properties of the tree might influence the error rates, we modify the branch lengths
# For this, we built four additional variants of each phylogenetic tree:
# Grafentree: a tree with branch lengths based on Grafen's method (all tips equidistant from root, branch length depends onnumber of nodes between root and tip)
# Onetree: a tree with all branch lengths set to have the same length of one
# Earlytree: a tree with early diversification and long branches leading to the tips
# Latetree: a tree with recent diversification and long branches between clades

source("https://raw.githubusercontent.com/dieterlukas/CulturalMacroevolution_Simulation/master/Code/SimulateCulturalEvolution_ModifyPhylogeny.R")


#------------------------------------------------------------------------------------------

# To determine whether properties of the sample might influence the error rates, we split the societies into different clades
# Societies in the largest lineage are classified as Clade A

source("https://raw.githubusercontent.com/dieterlukas/CulturalMacroevolution_Simulation/master/Code/SimulateCulturalEvolution_SubsetClade.R")


#------------------------------------------------------------------------------------------
# For the analyses to determine whether changes occurred differently in different parts of the tree, there are two specifications
# The first one takes the information from above to contrast changes in the largest clade (CladeA) from changes in the remaining part of the tree
# The second one takes information on an ecological variables to classify societies into two groups: 
# for the American societies this is whether they primarily use a forest habitat or not
# for the PamaNyungan societies this is whether they are hunter-gatherers or food-producers

# To determine whether changes occurred differently, we need to match this information from the societies to the tree
# For the clade-based distinction, this labels all branches within the clade one way an all other branches another way
# For the ecology-based distinction, branches are labelled based on a phylogenetic reconstruction of the most likely history of the ecological variable

source("https://raw.githubusercontent.com/dieterlukas/CulturalMacroevolution_Simulation/master/Code/SimulateCulturalEvolution_MatchPredictors.R")

#------------------------------------------------------------------------------------------
# This completes the preparation of the tree and data
#------------------------------------------------------------------------------------------


# We now select the correct settings for the specific run of the simulations

# The following selects the dataset, and sets up loops to repeat simulations 
# across all variants of the trees and using the different drift models


if(Option=="WNAI") {data<-WNAIdata}
if(Option=="PamaNyungan") {data<-PamaNyungandata}
if(Option=="WNAI") {locationdata<-WNAIlocations}
if(Option=="PamaNyungan") {locationdata<-PamaNyunganlocations}


# For the WNAI, the information to build the phylogeny does not provide branch lenghts. Therefore, we only
# have the phylogenies with the modified branch lengths

if(Option=="WNAI") {tree_variants <- c("Grafentree","Earlytree","Latetree")}
if(Option=="PamaNyungan") {tree_variants <- c("PamaNyungantree","Grafentree","Earlytree","Latetree")}

# We need a counter to keep track of the progression of the simulations and store the results
# in the correct line of the output file. 
counter<-1

# To make it easier to follow the process of the simulation we suppress warnings
oldw <- getOption("warn")
options(warn = -1)

# Some of the analyses use Bayesian regressions as in the package MCMCglmm, which require a prior
prior1 = list(R = list(V = 10, fix = 1), G = list(G1 = list(V = 1, nu = 0.002)))


# The simulations iterate across the different trees, across different subsets of the data, and with or without
# horizontal transmission. For each type of trait and inference, there are therefore multiple outputs.

# Change the number of repetitions to modify the number of indepedent simulations of each of these variants
# with a given drift model on a given tree are being performed and analysed
# This will influence how long the simulations take: with 10 repetitions, one set of analyses (e.g. discrete trait with no direction)
# can take ~6 hours on a standard desktop computer.
repetitions<-10


#-------------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------

# Each simulation loops over the three (WNAI) or four (PamaNyungan) different phylogenies. 

# For each phylogeny, setting and repetition, the following steps occur:
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


# For each sample (within each setting and repetition), analyses are performed to get the two types of error rates
# for the three different inferences. For each inference, one analysis model has an evolutionary model that assumes
# that there was no signal while another analysis model has an evolutionary model that assumes that there is no signal.
# The likelihoods of each model are saved. For the simulations with no signal, a false positive error occurs
# when the likelihood of the model with signal is significantly better than the likelihood of the model without a signal.
# For the simulations with a signal, a false negative error occurs when the likelihood of the model with signal
# is not significantly better than the likelihood of the model without a signal.

# There is an analysis script the shows how to potentially analyse the output from the simulations.



#-------------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------


# There are six different options for the simulations: traits are either discrete or continuous; and in each case
# changes might occur randomly, in a direction, or be influenced by where in the phylogeny they occur.

# Discrete traits reflect variables that can be split into two or more discrete categories.
# This reflects a large number of cultural traits, which are either classified as present/absent or
# into different separate categories. The simulated traits have either two or four variants.

# Continuous traits reflect variables that can be measured along a continuum. 
# This reflects traits such as body size or population density. The simulated traits start from a value of 0,
# and from there can increase or decrease.

# Drift code simulates changes in an arbitrary cultural trait across a known phylogeny, assuming that all 
# changes are equally likely, in all parts of the tree. This reflects a trait that changed according to drift.
# The resulting distribution of the different states of the trait are analysed with a variety of methods
# that are frequently used in biological and cultural phylogenetics. A false positive error occurs
# if such an analysis wrongly infers that one of the variants of the trait has been under selection or
# that changes in the trait occurred at different rates in different parts of the tree.

source("https://raw.githubusercontent.com/dieterlukas/CulturalMacroevolution_Simulation/master/Code/SimulateCulturalEvolution_Continuous_Drift.R")

source("https://raw.githubusercontent.com/dieterlukas/CulturalMacroevolution_Simulation/master/Code/SimulateCulturalEvolution_Discrete_Drift.R")


# Directional code simulates changes in an arbitrary cultural trait across a known phylogeny, assuming that some
# changes are more likely than others. For the discrete traits, this assumes that transitions to and from one state
# are more common than other transitions. For the continuous traits, this assumes that increases are more likely to
# occur than decreases. There are further options that could be changed directly in the respective code files.

source("https://raw.githubusercontent.com/dieterlukas/CulturalMacroevolution_Simulation/master/Code/SimulateCulturalEvolution_Continuous_Directional.R")

source("https://raw.githubusercontent.com/dieterlukas/CulturalMacroevolution_Simulation/master/Code/SimulateCulturalEvolution_Discrete_Directional.R")


# Clade difference code simulates changes in an arbitrary cultural trait across a known phylogeny, assuming that 
# certain changes are more likely in some parts of the tree than in others. Transitions are either different in 
# Clade A compared to the rest of the tree, or in lineages with a particular ecology compared to the rest of the tree.


source("https://raw.githubusercontent.com/dieterlukas/CulturalMacroevolution_Simulation/master/Code/SimulateCulturalEvolution_Continuous_Directional.R")

source("https://raw.githubusercontent.com/dieterlukas/CulturalMacroevolution_Simulation/master/Code/SimulateCulturalEvolution_Discrete_Directional.R")

















# Restore the usual patterns of warning in R
options(warn = oldw)


