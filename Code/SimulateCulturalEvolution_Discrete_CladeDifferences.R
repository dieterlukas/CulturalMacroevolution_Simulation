

#------------------------------------------------------------------------------------------
# Simulations for the manuscript - selection
#------------------------------------------------------------------------------------------

# These simulations and their outcomes are described in:
# Lukas, D., Towner, M., & Mulder, M. B. (2020). The Potential to Infer the Historical Pattern of Cultural Macroevolution as Illustrated by the Western North American Indian Societies.
# https://osf.io/preprints/socarxiv/tjvgy/




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





predictor_variants<-c("CladeA","Ecology")


#This will start the loop; 
#because of the large number of analyses for each simulation (different techniques, various subsamples), this takes noticeable computer time

for (tree_variant in 1:length(tree_variants)) {
  
  tree_used<-tree_variants[tree_variant]
  if(tree_used=="PamaNyungantree") currenttree<-PamaNyungantree
  if(tree_used=="Grafentree") currenttree<-Grafentree
  if(tree_used=="Onetree") currenttree<-Onetree
  if(tree_used=="Earlytree") currenttree<-Earlytree
  if(tree_used=="Latetree") currenttree<-Latetree

  for (predictor_variant in 1:length(predictor_variants)) {
    
    predictorvariant_used<-predictor_variants[predictor_variant]
    if(selectionvariant_used=="CladeA") current_predictorvariant<-simmapattemptCladeA
    if(selectionvariant_used=="Ecology") current_predictorvariant<-simmapattemptEcology
    
  for (selection_variant in 1:length(selection_variants)) {
    
    selectionvariant_used<-selection_variants[selection_variant]
    if(selectionvariant_used=="Q_onlyandfastInEcology1") current_selectionvariant<-Q_onlyandfastInEcology1
    if(selectionvariant_used=="Q_moreandslowInEcology1") current_selectionvariant<-Q_moreandslowInEcology1
    if(selectionvariant_used=="Q_onlyandslowInEcology1") current_selectionvariant<-Q_onlyandslowInEcology1
    if(selectionvariant_used=="Q_moreandfastInEcology1") current_selectionvariant<-Q_moreandfastInEcology1
        

    for (repetition in 1:repetitions) {
      
      
      #Simulate the data with the specified tree and the specified drift model
      
      simulatedtips_neutral<-sim.multiMk(tree=current_predictorvariant,Q=current_selectionvariant,anc="a")
      
      
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
          resultphylosiglambda<-phylosig(currenttree,numeric_simulatedtipds_neutral,method="lambda") 
          selection_results[counter,12]<-resultphylosiglambda$lambda
          resultphylosigK<-phylosig(currenttree,numeric_simulatedtipds_neutral,method="K") 
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
          resultphylosiglambda<-phylosig(currenttree,numeric_simulatedtipds_neutral,method="lambda") 
          selection_results[counter,12]<-resultphylosiglambda$lambda
          resultphylosigK<-phylosig(currenttree,numeric_simulatedtipds_neutral,method="K") 
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
          resultphylosiglambda<-phylosig(currenttree,numeric_simulatedtipds_neutral,method="lambda") 
          selection_results[counter,12]<-resultphylosiglambda$lambda
          resultphylosigK<-phylosig(currenttree,numeric_simulatedtipds_neutral,method="K") 
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
          resultphylosiglambda<-phylosig(currenttree,numeric_simulatedtipds_neutral,method="lambda") 
          selection_results[counter,12]<-resultphylosiglambda$lambda
          resultphylosigK<-phylosig(currenttree,numeric_simulatedtipds_neutral,method="K") 
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
        
        
        
        print(counter)
        print(repetition)
        print(tree_used)
        print(selection_variant)
        
        
      } #only need to run if there is more than one variant  
      
    } # end of the repetition loop
    
    if(Option=="WNAI") write.csv(drift_results[1:counter,],file="SimulateCulture_EcologyResults_Discrete_WNAI.csv")  
    if(Option=="PamaNyungan") write.csv(drift_results[1:counter,],file="SimulateCulture_EcologyResults_Discrete_PamaNyungan.csv")  
    
  }  # end of the drift variants loop
  
  } # end of the predictor variant loop  
  
} #end of the tree variant loop

