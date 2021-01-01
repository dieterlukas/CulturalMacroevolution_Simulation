

#------------------------------------------------------------------------------------------
# Simulations for the manuscript - continuous drift
#------------------------------------------------------------------------------------------

# These simulations and their outcomes are described in:
# Lukas, D., Towner, M., & Mulder, M. B. (2020). The Potential to Infer the Historical Pattern of Cultural Macroevolution as Illustrated by the Western North American Indian Societies.
# https://osf.io/preprints/socarxiv/tjvgy/



#------------------------------------------------------------------------------------------
# This set of simulations and analyses aims to assess the potential false positive error rate 
# of phylogenetic reconstructions of the historical evolution of cultural traits

# The code simulates changes in an arbitrary cultural trait across a known phylogeny, assuming that all 
# changes are equally likely, in all parts of the tree. This reflects a trait that changed according to drift.
# The resulting distribution of the different states of the trait are analysed with a variety of methods
# that are frequently used in biological and cultural phylogenetics. A false positive error occurs
# if such an analysis wrongly infers that one of the variants of the trait has been under selection or
# that changes in the trait occurred at different rates in different parts of the tree.

# The analyses examine three potential false positive errors:
# 1) wrongly identifying selection on the trait:
#         determines whether an analysis wrongly supports the inference that changes to or from one state
#         have occurred significantly more likely than other changes - even though the simulation is set up
#         such that changes between all states are equally likely
# 2) wrongly identifying lineage differences in changes in the trait:
#         determines whether an analysis wrongly supports the inference that changes between states
#         occurred more or less frequently in one part of the tree compared to the rest of the tree - even
#         though the simulation is set up such that changes occur at the same rate throughout the tree
# 3) wronlgy identifying ecological associations of the distribution of the trait:
#         determines whether an analysis wrongly supports that changes to or from one state have occurred
#         significantly more often in lineages that have a particular ecology - even though the simulation
#         is set up such that changes in the trait are not associated with ecological differences


# The setup includes further modifications to assess whether the false positive error rates are influenced by
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


#We prepare a data frame to store all the output
#Each row will have the results from a single simulation (with a given simulation model on a given tree)
#The columns contain basic information on the settings and the output from the various reconstruction methods
drift_results<-matrix(data=NA,nrow=500000,ncol=7)
colnames(drift_results)<-c("Tree","Repetition","Sample","pvalue_ecology","pvalue_clade","variance_trait","pvalue_trend")




#This will start the loop; 
#because of the large number of analyses for each simulation (different techniques, various subsamples), this takes noticeable computer time

for (tree_variant in 1:length(tree_variants)) {
  
  tree_used<-tree_variants[tree_variant]
  if(tree_used=="PamaNyungantree") currenttree<-PamaNyungantree
  if(tree_used=="Grafentree") currenttree<-Grafentree
  if(tree_used=="Onetree") currenttree<-Onetree
  if(tree_used=="Earlytree") currenttree<-Earlytree
  if(tree_used=="Latetree") currenttree<-Latetree
  
    
    for (repetition in 1:repetitions) {
      
      
      #Simulate the data with the specified tree and the specified drift model
      

      #From the OUwie.sim settings:
      # a. Single rate Brownian motion (BM1): alpha=c(1e-10,1e-10); sigma.sq=c(0.45,0.45); theta0=1.0; theta=c(0,0).
      
      # b. Brownian motion with different rate parameters for each state on a tree (BMS): alpha=c(1e-10,1e-10) sigma.sq=c(0.45,0.90); theta0=1.0; theta=c(0,0).
      
      # c. Ornstein Uhlenbeck with a single optimum for all species (OU1): alpha=c(0.1,0.1); sigma.sq=c(0.9,0.9); theta0=1; theta=c(1.0,1.0).
      
      # d. Ornstein Uhlenbeck model that assumes different state means and a single alpha and sigma^2 (OUM): alpha=c(1.0,1.0); sigma.sq=c(0.45,0.45); theta0=1.0; theta=c(1.0,2.0).
      
      # e. Ornstein Uhlenbeck model that assumes different state means and multiple sigma^2 (OUMV): alpha=c(1.0,1.0); sigma.sq=c(0.45,0.90); theta0=1.0; theta=c(1.0,2.0).
      
      # f. Ornstein Uhlenbeck model that assumes different state means and multiple alpha (OUMA): alpha=c(1.0,0.5); sigma.sq=c(0.45,0.45); theta0=1.0; theta=c(1.0,2.0).
      
      # g. Ornstein Uhlenbeck model that assumes different state means and multiple sigma^2 and alpha (OUMVA): alpha=c(1.0,0.5); sigma.sq=c(0.45,0.9); theta0=1.0; theta=c(1.0,2.0).
      
      
      
      simulatedtips_direction<-OUwie.sim(phy=simmapattemptEcology,simmap.tree = TRUE,alpha=c(1e-10,1e-10), sigma.sq=c(0.45,0.45), theta0=1, theta=c(0,0))
      
      simulatedtips_neutral<-simulatedtips_direction$X
      names(simulatedtips_neutral)<-simulatedtips_direction$Genus_species
      

      
      
      
      
      

      internalnodes<-fastAnc(currenttree,simulatedtips_neutral)
      
       
      #store the data in case something breaks
      simulatedtips_problem<-simulatedtips_neutral
      
      
      #Convert the data for the analyses
      #The matrix for the selection reconstruction assumes that the first value is the one under selection
      #We are making the assumption that if there is selection, it would have acted on the variant present among the largest number of societies
      #We therefore assign the value that appears most commonly among the simulated tip values as the first value
      
      
        
        joined_ecology<-data.frame(names(simulatedtips_neutral),Ecologymember,as.numeric(simulatedtips_neutral))
        colnames(joined_ecology)<-c("animal","subsistence","simulated")
        
        
        INphylo<-inverseA(currenttree)
        MCMC_regression_ecology<-MCMCglmm(simulated ~ subsistence,random=~animal,data=joined_ecology,pedigree=INphylo$pedigree,nitt=13000,prior=prior1,verbose=F,burnin=3000)
        
        
        joined_lineage<-data.frame(as.numeric(simulatedtips_neutral),cladeAmember,names(simulatedtips_neutral))
        colnames(joined_lineage)<-c("simulated","cladeA","animal")
        
        
        MCMC_regression_lineage<-MCMCglmm(simulated ~ cladeA,random=~animal,data=joined_lineage,pedigree=INphylo$pedigree,nitt=13000,prior=prior1,verbose=F,burnin=3000)
        
        
          
        
        #Start filling in the data in the respective column
        drift_results[counter,1]<-tree_used
        drift_results[counter,2]<-repetition
        
        drift_results[counter,3]<-"full"
        
        if(length(unique(joined_ecology$subsistence))==2) drift_results[counter,4]<-summary(MCMC_regression_ecology)$solutions[2,5]
        
        if(length(unique(joined_lineage$cladeA))==2) drift_results[counter,5]<-summary(MCMC_regression_lineage)$solutions[2,5]
        
        drift_results[counter,6]<-var(simulatedtips_neutral)
        
        
        
        reconstructionresults<-list() # setting up an empty list for the results
        mods<-c("BM1", "OU1") # We set up a vector containing the names of the models we want to run and use this in a simple loop.
        
        for(i in 1:length(mods)){
          reconstructionresults[[i]]<-OUwie(currenttree, joined_ecology, model=mods[i], simmap.tree=FALSE, root.station=TRUE,warn=FALSE,quiet = TRUE)
        }
        
        
        
        drift_results[counter,7]<-pchisq( (2*reconstructionresults[[1]]$loglik-reconstructionresults[[2]]$loglik)   , 1 )
       
        
        counter<-counter+1
        
        
        
#--#--#--#--#--#--#        #--#--#--#--#--#--#
        
        #Repeat the analyses with the different subsamples
        
        
           
        
        
#--#--#--#--#--#--#        #-#-#-#-#-#-#-#-#-#-#-#
        
        # Subset the tree and the data to assume that 25% of data is randomly missing
        
        
          
          ThreeQuarterClade <- currenttree$tip.label
          ThreeQuarterClade <- sample(ThreeQuarterClade)
          ThreeQuarterClade <- ThreeQuarterClade[1:129]
          
          ThreeQuartertree<-keep.tip(currenttree,ThreeQuarterClade)
          
          ThreeQuartersimulatedtips<-simulatedtips_neutral[names(simulatedtips_neutral) %in% ThreeQuarterClade]
          
        
          
          joined_ecology<-data.frame(names(ThreeQuartersimulatedtips),Ecologymember[names(Ecologymember) %in% names(ThreeQuartersimulatedtips)],as.numeric(ThreeQuartersimulatedtips))
          colnames(joined_ecology)<-c("animal","subsistence","simulated")
          
          
          INphylo<-inverseA(ThreeQuartertree)
          MCMC_regression_ecology<-MCMCglmm(simulated ~ subsistence,random=~animal,data=joined_ecology,pedigree=INphylo$pedigree,nitt=13000,prior=prior1,verbose=F,burnin=3000)
          
          
          
          
          joined_lineage<-data.frame(names(ThreeQuartersimulatedtips)  ,cladeAmember[names(cladeAmember) %in% names(ThreeQuartersimulatedtips)],as.numeric(ThreeQuartersimulatedtips))
          colnames(joined_lineage)<-c("animal","cladeA","simulated")
          
          
          MCMC_regression_lineage<-MCMCglmm(simulated ~ cladeA,random=~animal,data=joined_lineage,pedigree=INphylo$pedigree,nitt=13000,prior=prior1,verbose=F,burnin=3000)
          
          
          
          
          #Start filling in the data in the respective column
          drift_results[counter,1]<-tree_used
          drift_results[counter,2]<-repetition
          
          drift_results[counter,3]<-"threequarter"
          
          if(length(unique(joined_ecology$subsistence))==2) drift_results[counter,4]<-summary(MCMC_regression_ecology)$solutions[2,5]
          
          if(length(unique(joined_lineage$cladeA))==2) drift_results[counter,5]<-summary(MCMC_regression_lineage)$solutions[2,5]
          
          drift_results[counter,6]<-var(simulatedtips_neutral)
          
          
          reconstructionresults<-list() # setting up an empty list for the results
          mods<-c("BM1", "OU1") # We set up a vector containing the names of the models we want to run and use this in a simple loop.
          
          for(i in 1:length(mods)){
            reconstructionresults[[i]]<-OUwie(ThreeQuartertree, joined_ecology, model=mods[i], root.station=TRUE,warn=FALSE,quiet = TRUE)
          }
          
          
          
          drift_results[counter,7]<-pchisq( (2*reconstructionresults[[1]]$loglik-reconstructionresults[[2]]$loglik)   , 1 )
          
          
          counter<-counter+1
          
        
        
        
        
        
        
        #-#-#-#-#-#-#-#-#-#-#-#
        
        # Subset the tree and the data to assume that 50% of data is randomly missing
        
        
          
          HalfClade <- currenttree$tip.label
          HalfClade <- sample(HalfClade)
          HalfClade <- HalfClade[1:round(nrow(data)/2,0)]
          
          Halftree<-keep.tip(currenttree,HalfClade)
          
          Halfsimulatedtips<-simulatedtips_neutral[names(simulatedtips_neutral) %in% HalfClade]
          
          
          
          joined_ecology<-data.frame(names(Halfsimulatedtips) ,Ecologymember[names(Ecologymember) %in% names(Halfsimulatedtips)],as.numeric(Halfsimulatedtips))
          colnames(joined_ecology)<-c("animal","subsistence","simulated")
          
          
          INphylo<-inverseA(Halftree)
          MCMC_regression_ecology<-MCMCglmm(simulated ~ subsistence,random=~animal,data=joined_ecology,pedigree=INphylo$pedigree,nitt=13000,prior=prior1,verbose=F,burnin=3000)
          
          
          
          
          joined_lineage<-data.frame(names(Halfsimulatedtips)  ,cladeAmember[names(cladeAmember) %in% names(Halfsimulatedtips)],as.numeric(Halfsimulatedtips))
          colnames(joined_lineage)<-c("animal","cladeA","simulated")
          
          
          MCMC_regression_lineage<-MCMCglmm(simulated ~ cladeA,random=~animal,data=joined_lineage,pedigree=INphylo$pedigree,nitt=13000,prior=prior1,verbose=F,burnin=3000)
          
          
          
          
          #Start filling in the data in the respective column
          drift_results[counter,1]<-tree_used
          drift_results[counter,2]<-repetition
          
          drift_results[counter,3]<-"half"
          
          if(length(unique(joined_ecology$subsistence))==2) drift_results[counter,4]<-summary(MCMC_regression_ecology)$solutions[2,5]
          
          if(length(unique(joined_lineage$cladeA))==2) drift_results[counter,5]<-summary(MCMC_regression_lineage)$solutions[2,5]
          
          drift_results[counter,6]<-var(simulatedtips_neutral)
          
          reconstructionresults<-list() # setting up an empty list for the results
          mods<-c("BM1", "OU1") # We set up a vector containing the names of the models we want to run and use this in a simple loop.
          
          for(i in 1:length(mods)){
            reconstructionresults[[i]]<-OUwie(Halftree, joined_ecology, model=mods[i],  root.station=TRUE,warn=FALSE,quiet = TRUE)
          }
          
          
          
          drift_results[counter,7]<-pchisq( (2*reconstructionresults[[1]]$loglik-reconstructionresults[[2]]$loglik)   , 1 )
          
          
          counter<-counter+1
          
        
        
        
        
        
        #-#-#-#-#-#-#-#-#-#-#-#
        
        # State dependent loss of samples - 50% from positive values 
        
          
          randomdeletion<-function(x) (x*rbinom(1,1,0.5))
          
          
          
          PositiveValues<-names(simulatedtips_neutral[simulatedtips_neutral>0])
          NegativeValues<-names(simulatedtips_neutral[simulatedtips_neutral<0])
          PositiveValues <- sample(PositiveValues)
          PositiveValues <- PositiveValues[1:round(nrow(data)/2,0)]
          RetainedValues<-c(NegativeValues,PositiveValues)
          

            
            Rarelostsimulatedtips<-simulatedtips_neutral[names(simulatedtips_neutral) %in% RetainedValues]
            Rarelost<-names(Rarelostsimulatedtips)
            Rarelosttree<-keep.tip(currenttree,Rarelost)
          
          
          
          
            joined_ecology<-data.frame(names(Rarelostsimulatedtips) ,Ecologymember[names(Ecologymember) %in% names(Rarelostsimulatedtips)],as.numeric(Rarelostsimulatedtips))
            colnames(joined_ecology)<-c("animal","subsistence","simulated")
            
            
            INphylo<-inverseA(Rarelosttree)
            MCMC_regression_ecology<-MCMCglmm(simulated ~ subsistence,random=~animal,data=joined_ecology,pedigree=INphylo$pedigree,nitt=13000,prior=prior1,verbose=F,burnin=3000)
            
            
            
            joined_lineage<-data.frame(as.numeric(Rarelostsimulatedtips),cladeAmember[names(cladeAmember) %in% names(Rarelostsimulatedtips)],names(Rarelostsimulatedtips))
            colnames(joined_lineage)<-c("simulated","cladeA","animal")
            
            
            MCMC_regression_lineage<-MCMCglmm(simulated ~ cladeA,random=~animal,data=joined_lineage,pedigree=INphylo$pedigree,nitt=13000,prior=prior1,verbose=F,burnin=3000)
            
            
            
            
            #Start filling in the data in the respective column
            drift_results[counter,1]<-tree_used
            drift_results[counter,2]<-repetition
            
            drift_results[counter,3]<-"positiveloss"
            
            if(length(unique(joined_ecology$subsistence))==2) drift_results[counter,4]<-summary(MCMC_regression_ecology)$solutions[2,5]
            
            if(length(unique(joined_lineage$cladeA))==2) drift_results[counter,5]<-summary(MCMC_regression_lineage)$solutions[2,5]
            
            drift_results[counter,6]<-var(simulatedtips_neutral)
          
            reconstructionresults<-list() # setting up an empty list for the results
            mods<-c("BM1", "OU1") # We set up a vector containing the names of the models we want to run and use this in a simple loop.
            
            for(i in 1:length(mods)){
              reconstructionresults[[i]]<-OUwie(Rarelosttree, joined_ecology, model=mods[i], root.station=TRUE,warn=FALSE,quiet = TRUE)
            }
            
            
            
            drift_results[counter,7]<-pchisq( (2*reconstructionresults[[1]]$loglik-reconstructionresults[[2]]$loglik)   , 1 )
            
            
          counter<-counter+1
          
        
        
        
        
        
        
        
        
        # Loss according to ecoregion (50% of those living in forests or hunting/gathering)
        
        
        # Subset the tree and the data to assume that 50% of data is randomly missing
        
        
        
          
          randomdeletion<-function(x) (x*rbinom(1,1,0.5))
          Ecologytips<-sapply(as.numeric(Ecologymember),randomdeletion)
          
          
          if(Option=="WNAI") Ecologylosstips<-simulatedtips_neutral[Ecologytips==0]
          if(Option=="PamaNyungan") Ecologylosstips<-simulatedtips_neutral[Ecologytips==1]
          
          Ecologylosttree<-keep.tip(currenttree,names(Ecologylosstips))
          
          
          
          
          joined_ecology<-data.frame(names(Ecologylosstips)  ,Ecologymember[names(Ecologymember) %in% names(Ecologylosstips)],as.numeric(Ecologylosstips))
          colnames(joined_ecology)<-c("animal","subsistence","simulated")
          
          
          INphylo<-inverseA(Ecologylosttree)
          MCMC_regression_ecology<-MCMCglmm(simulated ~ subsistence,random=~animal,data=joined_ecology,pedigree=INphylo$pedigree,nitt=13000,prior=prior1,verbose=F,burnin=3000)
          
          
          
          
          joined_lineage<-data.frame(names(Ecologylosstips)  ,cladeAmember[names(cladeAmember) %in% names(Ecologylosstips)],as.numeric(Ecologylosstips))
          colnames(joined_lineage)<-c("animal","cladeA","simulated")
          
          
          MCMC_regression_lineage<-MCMCglmm(simulated ~ cladeA,random=~animal,data=joined_lineage,pedigree=INphylo$pedigree,nitt=13000,prior=prior1,verbose=F,burnin=3000)
          
          
          
          
          #Start filling in the data in the respective column
          drift_results[counter,1]<-tree_used
          drift_results[counter,2]<-repetition
          
          drift_results[counter,3]<-"ecologyloss"
          
          if(length(unique(joined_ecology$subsistence))==2) drift_results[counter,4]<-summary(MCMC_regression_ecology)$solutions[2,5]
          
          if(length(unique(joined_lineage$cladeA))==2) drift_results[counter,5]<-summary(MCMC_regression_lineage)$solutions[2,5]
          
          drift_results[counter,6]<-var(simulatedtips_neutral)
          
          reconstructionresults<-list() # setting up an empty list for the results
          mods<-c("BM1", "OU1") # We set up a vector containing the names of the models we want to run and use this in a simple loop.
          
          if(sum(joined_ecology[,2])==nrow(joined_ecology)){joined_ecology[1,2]<-0}
          
          for(i in 1:length(mods)){
            reconstructionresults[[i]]<-OUwie(Ecologylosttree, joined_ecology, model=mods[i], root.station=TRUE,warn=FALSE,quiet = TRUE)
          }
          
          
          
          drift_results[counter,7]<-pchisq( (2*reconstructionresults[[1]]$loglik-reconstructionresults[[2]]$loglik)   , 1 )
          
          
          counter<-counter+1
          
        
        
        

        
        
                
#--#--#--#--#--#--#        #-#-#-#-#-#-#-#-#-#-#-#-#
        
        
        # Horizontal transfer, 10% of societies change, randomly adopt another state
        
          
          
          societiestochange<-sample(simulatedtips_neutral,round(0.1*nrow(data),0))
          societiestocopy<-sample(simulatedtips_neutral,round(0.1*nrow(data),0))
          
          
          horizontal10tips<-simulatedtips_neutral
          
          for (i in 1:round(0.1*nrow(data),0)) {
          currentchange<-  societiestochange[i]
          horizontal10tips[names(horizontal10tips) %in% names(currentchange)]<-societiestocopy[i]
          }
        
          
          joined_ecology<-data.frame(names(horizontal10tips) ,Ecologymember,as.numeric(horizontal10tips))
          colnames(joined_ecology)<-c("animal","subsistence","simulated")
          
          
          INphylo<-inverseA(currenttree)
          MCMC_regression_ecology<-MCMCglmm(simulated ~ subsistence,random=~animal,data=joined_ecology,pedigree=INphylo$pedigree,nitt=13000,prior=prior1,verbose=F,burnin=3000)
          
          
          
          
          joined_lineage<-data.frame(as.numeric(horizontal10tips),cladeAmember,names(horizontal10tips))
          colnames(joined_lineage)<-c("simulated","cladeA","animal")
          
          
          MCMC_regression_lineage<-MCMCglmm(simulated ~ cladeA,random=~animal,data=joined_lineage,pedigree=INphylo$pedigree,nitt=13000,prior=prior1,verbose=F,burnin=3000)
          
          
          
          #Start filling in the data in the respective column
          drift_results[counter,1]<-tree_used
          drift_results[counter,2]<-repetition
          
          drift_results[counter,3]<-"horizontal10"
          
          if(length(unique(joined_ecology$subsistence))==2) drift_results[counter,4]<-summary(MCMC_regression_ecology)$solutions[2,5]
          
          if(length(unique(joined_lineage$cladeA))==2) drift_results[counter,5]<-summary(MCMC_regression_lineage)$solutions[2,5]
          
          drift_results[counter,6]<-var(simulatedtips_neutral)
         
          reconstructionresults<-list() # setting up an empty list for the results
          mods<-c("BM1", "OU1") # We set up a vector containing the names of the models we want to run and use this in a simple loop.
          
          for(i in 1:length(mods)){
            reconstructionresults[[i]]<-OUwie(currenttree, joined_ecology, model=mods[i], simmap.tree=FALSE, root.station=TRUE,warn=FALSE,quiet = TRUE)
          }
          
          
          
          drift_results[counter,7]<-pchisq( (2*reconstructionresults[[1]]$loglik-reconstructionresults[[2]]$loglik)   , 1 )
          
          
          counter<-counter+1
          

        
        
        
        
        
        # Horizontal transfer, 10% of societies change, adopt state from within their own clade
        
        
          
          societiestochange<-sample(simulatedtips_neutral[names(simulatedtips_neutral) %in% names(cladeAmember[cladeAmember==1])],round(0.1*nrow(data),0))
          societiestocopy<-sample(simulatedtips_neutral[names(simulatedtips_neutral) %in% names(cladeAmember[cladeAmember==1])],round(0.1*nrow(data),0))
          
          
          horizontal10tips<-simulatedtips_neutral
          
          for (i in 1:round(0.1*nrow(data),0)) {
            currentchange<-  societiestochange[i]
            horizontal10tips[names(horizontal10tips) %in% names(currentchange)]<-societiestocopy[i]
          }
          
          
          joined_ecology<-data.frame(names(horizontal10tips) ,Ecologymember,as.numeric(horizontal10tips))
          colnames(joined_ecology)<-c("animal","subsistence","simulated")
          
          
          INphylo<-inverseA(currenttree)
          MCMC_regression_ecology<-MCMCglmm(simulated ~ subsistence,random=~animal,data=joined_ecology,pedigree=INphylo$pedigree,nitt=13000,prior=prior1,verbose=F,burnin=3000)
          
          
          
          
          joined_lineage<-data.frame(as.numeric(horizontal10tips),cladeAmember,names(horizontal10tips))
          colnames(joined_lineage)<-c("simulated","cladeA","animal")
          
          
          MCMC_regression_lineage<-MCMCglmm(simulated ~ cladeA,random=~animal,data=joined_lineage,pedigree=INphylo$pedigree,nitt=13000,prior=prior1,verbose=F,burnin=3000)
          
          
          
          #Start filling in the data in the respective column
          drift_results[counter,1]<-tree_used
          drift_results[counter,2]<-repetition
          
          drift_results[counter,3]<-"horizontal10clade"
          
          if(length(unique(joined_ecology$subsistence))==2) drift_results[counter,4]<-summary(MCMC_regression_ecology)$solutions[2,5]
          
          if(length(unique(joined_lineage$cladeA))==2) drift_results[counter,5]<-summary(MCMC_regression_lineage)$solutions[2,5]
          
          drift_results[counter,6]<-var(simulatedtips_neutral)
          
          reconstructionresults<-list() # setting up an empty list for the results
          mods<-c("BM1", "OU1") # We set up a vector containing the names of the models we want to run and use this in a simple loop.
          
          for(i in 1:length(mods)){
            reconstructionresults[[i]]<-OUwie(currenttree, joined_ecology, model=mods[i], simmap.tree=FALSE, root.station=TRUE,warn=FALSE,quiet = TRUE)
          }
          
          
          
          drift_results[counter,7]<-pchisq( (2*reconstructionresults[[1]]$loglik-reconstructionresults[[2]]$loglik)   , 1 )
          
          
          counter<-counter+1
          

        
        
        # Horizontal transfer, 10% of societies change, adopt state from another society in same ecoregion
        
        
          
          societiestochange<-sample(simulatedtips_neutral[names(simulatedtips_neutral) %in% names(Ecologymember[Ecologymember==1])],round(0.1*nrow(data),0))
          societiestocopy<-sample(simulatedtips_neutral[names(simulatedtips_neutral) %in% names(Ecologymember[Ecologymember==1])],round(0.1*nrow(data),0))
          
          
          horizontal10tips<-simulatedtips_neutral
          
          for (i in 1:round(0.1*nrow(data),0)) {
            currentchange<-  societiestochange[i]
            horizontal10tips[names(horizontal10tips) %in% names(currentchange)]<-societiestocopy[i]
          }
          
          
          
          joined_ecology<-data.frame(names(horizontal10tips) ,Ecologymember,as.numeric(horizontal10tips))
          colnames(joined_ecology)<-c("animal","subsistence","simulated")
          
          
          INphylo<-inverseA(currenttree)
          MCMC_regression_ecology<-MCMCglmm(simulated ~ subsistence,random=~animal,data=joined_ecology,pedigree=INphylo$pedigree,nitt=13000,prior=prior1,verbose=F,burnin=3000)
          
          
          
          
          joined_lineage<-data.frame(as.numeric(horizontal10tips),cladeAmember,names(horizontal10tips))
          colnames(joined_lineage)<-c("simulated","cladeA","animal")
          
          
          MCMC_regression_lineage<-MCMCglmm(simulated ~ cladeA,random=~animal,data=joined_lineage,pedigree=INphylo$pedigree,nitt=13000,prior=prior1,verbose=F,burnin=3000)
          
          
          
          
          #Start filling in the data in the respective column
          drift_results[counter,1]<-tree_used
          drift_results[counter,2]<-repetition
          
          drift_results[counter,3]<-"horizontal10ecology"
          
          if(length(unique(joined_ecology$subsistence))==2) drift_results[counter,4]<-summary(MCMC_regression_ecology)$solutions[2,5]
          
          if(length(unique(joined_lineage$cladeA))==2) drift_results[counter,5]<-summary(MCMC_regression_lineage)$solutions[2,5]
          
          drift_results[counter,6]<-var(simulatedtips_neutral)
          
          reconstructionresults<-list() # setting up an empty list for the results
          mods<-c("BM1", "OU1") # We set up a vector containing the names of the models we want to run and use this in a simple loop.
          
          for(i in 1:length(mods)){
            reconstructionresults[[i]]<-OUwie(currenttree, joined_ecology, model=mods[i], simmap.tree=FALSE, root.station=TRUE,warn=FALSE,quiet = TRUE)
          }
          
          
          
          drift_results[counter,7]<-pchisq( (2*reconstructionresults[[1]]$loglik-reconstructionresults[[2]]$loglik)   , 1 )
          
          
          counter<-counter+1
          

          
          
          
          
          
          # Horizontal transfer, changes on 1% of branches
          
          nodeswith_ht<-sample(names(internalnodes),round(nrow(data)/100,0))
          nodeswith_ht<-sort(nodeswith_ht)
          
          for (nodestochange in 1:round(nrow(data)/100,0)) {
          currentnode<-nodeswith_ht[nodestochange]
          currenttips <- tips(currenttree, currentnode)
          ht_tree<-keep.tip(currenttree,currenttips)  
  
          after_ht_tips<-fastBM(tree=ht_tree,a=sample(internalnodes[names(internalnodes) != currentnode],1),sig2=1)
         
          ht_tips<-simulatedtips_neutral
          ht_tips[names(ht_tips) %in% names(after_ht_tips)]<-after_ht_tips
          }
   

          
          joined_ecology<-data.frame(names(ht_tips) ,Ecologymember,as.numeric(ht_tips))
          colnames(joined_ecology)<-c("animal","subsistence","simulated")
          
          
          INphylo<-inverseA(currenttree)
          MCMC_regression_ecology<-MCMCglmm(simulated ~ subsistence,random=~animal,data=joined_ecology,pedigree=INphylo$pedigree,nitt=13000,prior=prior1,verbose=F,burnin=3000)
          
          
          
          
          joined_lineage<-data.frame(as.numeric(ht_tips),cladeAmember,names(ht_tips))
          colnames(joined_lineage)<-c("simulated","cladeA","animal")
          
          
          MCMC_regression_lineage<-MCMCglmm(simulated ~ cladeA,random=~animal,data=joined_lineage,pedigree=INphylo$pedigree,nitt=13000,prior=prior1,verbose=F,burnin=3000)
          
          
          
          
          #Start filling in the data in the respective column
          drift_results[counter,1]<-tree_used
          drift_results[counter,2]<-repetition
          
          drift_results[counter,3]<-"horizontalnodes"
          
          if(length(unique(joined_ecology$subsistence))==2) drift_results[counter,4]<-summary(MCMC_regression_ecology)$solutions[2,5]
          
          if(length(unique(joined_lineage$cladeA))==2) drift_results[counter,5]<-summary(MCMC_regression_lineage)$solutions[2,5]
          
          drift_results[counter,6]<-var(simulatedtips_neutral)
          
          reconstructionresults<-list() # setting up an empty list for the results
          mods<-c("BM1", "OU1") # We set up a vector containing the names of the models we want to run and use this in a simple loop.
          
          for(i in 1:length(mods)){
            reconstructionresults[[i]]<-OUwie(currenttree, joined_ecology, model=mods[i], simmap.tree=FALSE, root.station=TRUE,warn=FALSE,quiet = TRUE)
          }
          
          
          
          drift_results[counter,7]<-pchisq( (2*reconstructionresults[[1]]$loglik-reconstructionresults[[2]]$loglik)   , 1 )
          
          
          counter<-counter+1  
          
          
        
          # Horizontal transfer, changes on 1% of branches in dominant ecology
          
          Ecologytree<-keep.tip(currenttree,names(Ecologymember[Ecologymember ==1]))
          
          Ecologyreconstructed<-fastAnc(currenttree,Ecologymember)
          Ecologyreconstructed<-round(Ecologyreconstructed,0)
          
          
          nodeswith_ht<-sample(names(internalnodes[Ecologyreconstructed==1]),round(nrow(data)/100,0))
          nodeswith_ht<-sort(nodeswith_ht)
          
          
          for (nodestochange in 1:round(nrow(data)/100,0)) {
            currentnode<-nodeswith_ht[nodestochange]
            currenttips <- tips(currenttree, currentnode)
            ht_tree<-keep.tip(currenttree,currenttips)  
            
            after_ht_tips<-fastBM(tree=ht_tree,a=sample(internalnodes[Ecologyreconstructed==1],1),sig2=1)
            
            ht_tips<-simulatedtips_neutral
            ht_tips[names(ht_tips) %in% names(after_ht_tips)]<-after_ht_tips
          }
          
          
          
          joined_ecology<-data.frame(names(ht_tips)  ,Ecologymember,as.numeric(ht_tips)) 
          colnames(joined_ecology)<-c("animal","subsistence","simulated")
          
          
          INphylo<-inverseA(currenttree)
          MCMC_regression_ecology<-MCMCglmm(simulated ~ subsistence,random=~animal,data=joined_ecology,pedigree=INphylo$pedigree,nitt=13000,prior=prior1,verbose=F,burnin=3000)
          
          
          
          joined_lineage<-data.frame(as.numeric(ht_tips),cladeAmember,names(ht_tips))
          colnames(joined_lineage)<-c("simulated","cladeA","animal")
          
          
          MCMC_regression_lineage<-MCMCglmm(simulated ~ cladeA,random=~animal,data=joined_lineage,pedigree=INphylo$pedigree,nitt=13000,prior=prior1,verbose=F,burnin=3000)
          
          
          
          
          #Start filling in the data in the respective column
          drift_results[counter,1]<-tree_used
          drift_results[counter,2]<-repetition
          
          drift_results[counter,3]<-"horizontalnodesecology"
          
          if(length(unique(joined_ecology$subsistence))==2) drift_results[counter,4]<-summary(MCMC_regression_ecology)$solutions[2,5]
          
          if(length(unique(joined_lineage$cladeA))==2) drift_results[counter,5]<-summary(MCMC_regression_lineage)$solutions[2,5]
          
          drift_results[counter,6]<-var(simulatedtips_neutral)
          
          reconstructionresults<-list() # setting up an empty list for the results
          mods<-c("BM1", "OU1") # We set up a vector containing the names of the models we want to run and use this in a simple loop.
          
          for(i in 1:length(mods)){
            reconstructionresults[[i]]<-OUwie(currenttree, joined_ecology, model=mods[i], simmap.tree=FALSE, root.station=TRUE,warn=FALSE,quiet = TRUE)
          }
          
          
          
          drift_results[counter,7]<-pchisq( (2*reconstructionresults[[1]]$loglik-reconstructionresults[[2]]$loglik)   , 1 )
          
          
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
          
          
          joined_ecology<-data.frame(names(horizontalNeighbor10tips) ,Ecologymember,as.numeric(horizontalNeighbor10tips))
          colnames(joined_ecology)<-c("animal","subsistence","simulated")
          
          
          INphylo<-inverseA(currenttree)
          MCMC_regression_ecology<-MCMCglmm(simulated ~ subsistence,random=~animal,data=joined_ecology,pedigree=INphylo$pedigree,nitt=13000,prior=prior1,verbose=F,burnin=3000)
          
          
          
          joined_lineage<-data.frame(as.numeric(horizontalNeighbor10tips),cladeAmember,names(horizontalNeighbor10tips))
          colnames(joined_lineage)<-c("simulated","cladeA","animal")
          
          
          MCMC_regression_lineage<-MCMCglmm(simulated ~ cladeA,random=~animal,data=joined_lineage,pedigree=INphylo$pedigree,nitt=13000,prior=prior1,verbose=F,burnin=3000)
          
          
          
          
          #Start filling in the data in the respective column
          drift_results[counter,1]<-tree_used
          drift_results[counter,2]<-repetition
          
          drift_results[counter,3]<-"horizontalneighbor"
          
          if(length(unique(joined_ecology$subsistence))==2) drift_results[counter,4]<-summary(MCMC_regression_ecology)$solutions[2,5]
          
          if(length(unique(joined_lineage$cladeA))==2) drift_results[counter,5]<-summary(MCMC_regression_lineage)$solutions[2,5]
          
          drift_results[counter,6]<-var(simulatedtips_neutral)
          
          reconstructionresults<-list() # setting up an empty list for the results
          mods<-c("BM1", "OU1") # We set up a vector containing the names of the models we want to run and use this in a simple loop.
          
          for(i in 1:length(mods)){
            reconstructionresults[[i]]<-OUwie(simmapattemptEcology, joined_ecology, model=mods[i], simmap.tree=TRUE, root.station=TRUE,warn=FALSE,quiet = TRUE)
          }
          
          
          
          drift_results[counter,7]<-pchisq( (2*reconstructionresults[[1]]$loglik-reconstructionresults[[2]]$loglik)   , 1 )
          
          
          counter<-counter+1
          
          
          
        print("runningcounter")
        print(counter)
        print("repetition")
        print(repetition)
        print("currenttree")
        print(tree_used)
        
        
       
      
    } # end of the repetition loop
    
  if(Option=="WNAI") write.csv(drift_results[1:counter,],file="SimulateCulture_NoDirectionalResults_Continuous_WNAI.csv")  
  if(Option=="PamaNyungan") write.csv(drift_results[1:counter,],file="SimulateCulture_NoDirectionalResults_Continuous_PamaNyungan.csv")  
  
} #end of the tree variant loop

