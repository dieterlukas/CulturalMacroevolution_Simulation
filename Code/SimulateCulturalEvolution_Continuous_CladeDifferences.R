

#------------------------------------------------------------------------------------------
# Simulations for the manuscript - continuous clade differences
#------------------------------------------------------------------------------------------

# These simulations and their outcomes are described in:
# Lukas, D., Towner, M., & Mulder, M. B. (2020). The Potential to Infer the Historical Pattern of Cultural Macroevolution as Illustrated by the Western North American Indian Societies.
# https://osf.io/preprints/socarxiv/tjvgy/



#We prepare a data frame to store all the output
#Each row will have the results from a single simulation (with a given simulation model on a given tree)
#The columns contain basic information on the settings and the output from the various reconstruction methods
drift_results<-matrix(data=NA,nrow=500000,ncol=14)
colnames(drift_results)<-c("Tree","Repetition","Sample","pvalue_ecology","pvalue_alternative","pvalue_clade","variance_trait","median_Ecology1","median_Ecology2","median_Alternative1","median_Alternative2","median_CladeA","median_OtherClades","predictor_variable")


predictor_variants<-c("simmapattemptCladeA","simmapattemptEcology","simmapattemptAlternative")

#This will start the loop; 
#because of the large number of analyses for each simulation (different techniques, various subsamples), this takes noticeable computer time

for (predictor in 1:length(predictor_variants)) {

predictor_variant_used<-predictor_variants[predictor]

for (tree_variant in 1:length(tree_variants)) {
  
  tree_used<-tree_variants[tree_variant]
  if(tree_used=="PamaNyungantree") currenttree<-PamaNyungantree
  if(tree_used=="Grafentree") currenttree<-Grafentree
  if(tree_used=="Onetree") currenttree<-Onetree
  if(tree_used=="Earlytree") currenttree<-Earlytree
  if(tree_used=="Latetree") currenttree<-Latetree
  
    
    for (repetition in 1:repetitions) {
      
      
      #Simulate the data with the specified tree and the specified drift model
      
      
      sig2<-setNames(c(2,2),c("0","1"))
      alpha<-setNames(c(4,4),c("0","1"))
      theta<-setNames(c(10,-10),c("0","1"))
      
      
      if(predictor_variant_used=="simmapattemptCladeA") simulatedtips_neutral<-multiOU(tree=simmapattemptCladeA,alpha,sig2,theta,internal=TRUE)
      if(predictor_variant_used=="simmapattemptEcology") simulatedtips_neutral<-multiOU(tree=simmapattemptEcology,alpha,sig2,theta,internal=TRUE)
      if(predictor_variant_used=="simmapattemptAlternative") simulatedtips_neutral<-multiOU(tree=simmapattemptAlternative,alpha,sig2,theta,internal=TRUE)
      
      simulatedtips_neutral<-simulatedtips_neutral[names(simulatedtips_neutral) %in% currenttree$tip.label]
      
      internalnodes<-fastAnc(currenttree,simulatedtips_neutral)
      
      
      #store the data in case something breaks
      simulatedtips_problem<-simulatedtips_neutral
      
      
      #Convert the data for the analyses
      #The matrix for the selection reconstruction assumes that the first value is the one under selection
      #We are making the assumption that if there is selection, it would have acted on the variant present among the largest number of societies
      #We therefore assign the value that appears most commonly among the simulated tip values as the first value
      
      
        

        
        
        
        joined_ecology<-data.frame(as.numeric(simulatedtips_neutral),Ecologymember,names(simulatedtips_neutral))
        colnames(joined_ecology)<-c("simulated","subsistence","animal")
        
        
        INphylo<-inverseA(currenttree)
        MCMC_regression_ecology<-MCMCglmm(simulated ~ subsistence,random=~animal,data=joined_ecology,pedigree=INphylo$pedigree,nitt=13000,prior=prior1,verbose=F,burnin=3000)
        
        
        joined_alternative<-data.frame(as.numeric(simulatedtips_neutral),Alternativemember,names(simulatedtips_neutral))
        colnames(joined_alternative)<-c("simulated","alternative","animal")
        
        
        INphylo<-inverseA(currenttree)
        MCMC_regression_alternative<-MCMCglmm(simulated ~ alternative,random=~animal,data=joined_alternative,pedigree=INphylo$pedigree,nitt=13000,prior=prior1,verbose=F,burnin=3000)
        
        
        
        joined_lineage<-data.frame(as.numeric(simulatedtips_neutral),cladeAmember,names(simulatedtips_neutral))
        colnames(joined_lineage)<-c("simulated","cladeA","animal")
        
        
        MCMC_regression_lineage<-MCMCglmm(simulated ~ cladeA,random=~animal,data=joined_lineage,pedigree=INphylo$pedigree,nitt=13000,prior=prior1,verbose=F,burnin=3000)
        
        
          
        
        #Start filling in the data in the respective column
        drift_results[counter,1]<-tree_used
        drift_results[counter,2]<-repetition
        
        drift_results[counter,3]<-"full"
        
        if(length(unique(joined_ecology$subsistence))==2) drift_results[counter,4]<-summary(MCMC_regression_ecology)$solutions[2,5]
        
        if(length(unique(joined_alternative$alternative))==2) drift_results[counter,5]<-summary(MCMC_regression_alternative)$solutions[2,5]
        
        if(length(unique(joined_lineage$cladeA))==2) drift_results[counter,6]<-summary(MCMC_regression_lineage)$solutions[2,5]
        
        drift_results[counter,7]<-var(simulatedtips_neutral)
        
        drift_results[counter,8]<-median(joined_ecology[joined_ecology$subsistence==0,]$simulated)
        drift_results[counter,9]<-median(joined_ecology[joined_ecology$subsistence==1,]$simulated)
        
        drift_results[counter,10]<-median(joined_alternative[joined_alternative$alternative==0,]$simulated)
        drift_results[counter,11]<-median(joined_alternative[joined_alternative$alternative==1,]$simulated)
        
        drift_results[counter,12]<-median(joined_lineage[joined_lineage$cladeA==0,]$simulated)
        drift_results[counter,13]<-median(joined_lineage[joined_lineage$cladeA==1,]$simulated)
        
        drift_results[counter,14]<-predictor_variant_used
        
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
          
        
          
          joined_ecology<-data.frame(as.numeric(ThreeQuartersimulatedtips),Ecologymember[names(Ecologymember) %in% names(ThreeQuartersimulatedtips)],names(ThreeQuartersimulatedtips))
          colnames(joined_ecology)<-c("simulated","subsistence","animal")
          
          
          INphylo<-inverseA(ThreeQuartertree)
          MCMC_regression_ecology<-MCMCglmm(simulated ~ subsistence,random=~animal,data=joined_ecology,pedigree=INphylo$pedigree,nitt=13000,prior=prior1,verbose=F,burnin=3000)
          
          
          joined_alternative<-data.frame(as.numeric(ThreeQuartersimulatedtips),Alternativemember[names(Alternativemember) %in% names(ThreeQuartersimulatedtips)],names(ThreeQuartersimulatedtips))
          colnames(joined_alternative)<-c("simulated","alternative","animal")
          
          
          INphylo<-inverseA(ThreeQuartertree)
          MCMC_regression_ecology<-MCMCglmm(simulated ~ alternative,random=~animal,data=joined_alternative,pedigree=INphylo$pedigree,nitt=13000,prior=prior1,verbose=F,burnin=3000)
          
          
          
          
          joined_lineage<-data.frame(as.numeric(ThreeQuartersimulatedtips),cladeAmember[names(cladeAmember) %in% names(ThreeQuartersimulatedtips)],names(ThreeQuartersimulatedtips))
          colnames(joined_lineage)<-c("simulated","cladeA","animal")
          
          
          MCMC_regression_lineage<-MCMCglmm(simulated ~ cladeA,random=~animal,data=joined_lineage,pedigree=INphylo$pedigree,nitt=13000,prior=prior1,verbose=F,burnin=3000)
          
          
          
          
          #Start filling in the data in the respective column
          drift_results[counter,1]<-tree_used
          drift_results[counter,2]<-repetition
          
          drift_results[counter,3]<-"quarterloss"
          
          if(length(unique(joined_ecology$subsistence))==2) drift_results[counter,4]<-summary(MCMC_regression_ecology)$solutions[2,5]
          
          if(length(unique(joined_alternative$alternative))==2) drift_results[counter,5]<-summary(MCMC_regression_alternative)$solutions[2,5]
          
          if(length(unique(joined_lineage$cladeA))==2) drift_results[counter,6]<-summary(MCMC_regression_lineage)$solutions[2,5]
          
          drift_results[counter,7]<-var(simulatedtips_neutral)
          
          drift_results[counter,8]<-median(joined_ecology[joined_ecology$subsistence==0,]$simulated)
          drift_results[counter,9]<-median(joined_ecology[joined_ecology$subsistence==1,]$simulated)
          
          drift_results[counter,10]<-median(joined_alternative[joined_alternative$alternative==0,]$simulated)
          drift_results[counter,11]<-median(joined_alternative[joined_alternative$alternative==1,]$simulated)
          
          drift_results[counter,12]<-median(joined_lineage[joined_lineage$cladeA==0,]$simulated)
          drift_results[counter,13]<-median(joined_lineage[joined_lineage$cladeA==1,]$simulated)
          
          drift_results[counter,14]<-predictor_variant_used
          
          counter<-counter+1
          
        
        
        
        
        
        
        #-#-#-#-#-#-#-#-#-#-#-#
        
        # Subset the tree and the data to assume that 50% of data is randomly missing
        
        
          
          HalfClade <- currenttree$tip.label
          HalfClade <- sample(HalfClade)
          HalfClade <- HalfClade[1:round(nrow(data)/2,0)]
          
          Halftree<-keep.tip(currenttree,HalfClade)
          
          Halfsimulatedtips<-simulatedtips_neutral[names(simulatedtips_neutral) %in% HalfClade]
          
          
          
          joined_ecology<-data.frame(as.numeric(Halfsimulatedtips),Ecologymember[names(Ecologymember) %in% names(Halfsimulatedtips)],names(Halfsimulatedtips))
          colnames(joined_ecology)<-c("simulated","subsistence","animal")
          
          
          INphylo<-inverseA(Halftree)
          MCMC_regression_ecology<-MCMCglmm(simulated ~ subsistence,random=~animal,data=joined_ecology,pedigree=INphylo$pedigree,nitt=13000,prior=prior1,verbose=F,burnin=3000)
          
          joined_alternative<-data.frame(as.numeric(Halfsimulatedtips),Alternativemember[names(Alternativemember) %in% names(Halfsimulatedtips)],names(Halfsimulatedtips))
          colnames(joined_alternative)<-c("simulated","alternative","animal")
          
          
          INphylo<-inverseA(Halftree)
          MCMC_regression_ecology<-MCMCglmm(simulated ~ alternative,random=~animal,data=joined_alternative,pedigree=INphylo$pedigree,nitt=13000,prior=prior1,verbose=F,burnin=3000)
          
          
          
          joined_lineage<-data.frame(as.numeric(Halfsimulatedtips),cladeAmember[names(cladeAmember) %in% names(Halfsimulatedtips)],names(Halfsimulatedtips))
          colnames(joined_lineage)<-c("simulated","cladeA","animal")
          
          
          MCMC_regression_lineage<-MCMCglmm(simulated ~ cladeA,random=~animal,data=joined_lineage,pedigree=INphylo$pedigree,nitt=13000,prior=prior1,verbose=F,burnin=3000)
          
          
          
          
          #Start filling in the data in the respective column
          drift_results[counter,1]<-tree_used
          drift_results[counter,2]<-repetition
          
          drift_results[counter,3]<-"halfloss"
          
          if(length(unique(joined_ecology$subsistence))==2) drift_results[counter,4]<-summary(MCMC_regression_ecology)$solutions[2,5]
          
          if(length(unique(joined_alternative$alternative))==2) drift_results[counter,5]<-summary(MCMC_regression_alternative)$solutions[2,5]
          
          if(length(unique(joined_lineage$cladeA))==2) drift_results[counter,6]<-summary(MCMC_regression_lineage)$solutions[2,5]
          
          drift_results[counter,7]<-var(simulatedtips_neutral)
          
          drift_results[counter,8]<-median(joined_ecology[joined_ecology$subsistence==0,]$simulated)
          drift_results[counter,9]<-median(joined_ecology[joined_ecology$subsistence==1,]$simulated)
          
          drift_results[counter,10]<-median(joined_alternative[joined_alternative$alternative==0,]$simulated)
          drift_results[counter,11]<-median(joined_alternative[joined_alternative$alternative==1,]$simulated)
          
          drift_results[counter,12]<-median(joined_lineage[joined_lineage$cladeA==0,]$simulated)
          drift_results[counter,13]<-median(joined_lineage[joined_lineage$cladeA==1,]$simulated)
          
          drift_results[counter,14]<-predictor_variant_used
          
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
          
          
          
          
            joined_ecology<-data.frame(as.numeric(Rarelostsimulatedtips),Ecologymember[names(Ecologymember) %in% names(Rarelostsimulatedtips)],names(Rarelostsimulatedtips))
            colnames(joined_ecology)<-c("simulated","subsistence","animal")
            
            
            INphylo<-inverseA(Rarelosttree)
            MCMC_regression_ecology<-MCMCglmm(simulated ~ subsistence,random=~animal,data=joined_ecology,pedigree=INphylo$pedigree,nitt=13000,prior=prior1,verbose=F,burnin=3000)
            
            joined_alternative<-data.frame(as.numeric(Rarelostsimulatedtips),Alternativemember[names(Alternativemember) %in% names(Rarelostsimulatedtips)],names(Rarelostsimulatedtips))
            colnames(joined_alternative)<-c("simulated","alternative","animal")
            
            
            INphylo<-inverseA(Rarelosttree)
            MCMC_regression_ecology<-MCMCglmm(simulated ~ alternative,random=~animal,data=joined_alternative,pedigree=INphylo$pedigree,nitt=13000,prior=prior1,verbose=F,burnin=3000)
            
            
            joined_lineage<-data.frame(as.numeric(Rarelostsimulatedtips),cladeAmember[names(cladeAmember) %in% names(Rarelostsimulatedtips)],names(Rarelostsimulatedtips))
            colnames(joined_lineage)<-c("simulated","cladeA","animal")
            
            
            MCMC_regression_lineage<-MCMCglmm(simulated ~ cladeA,random=~animal,data=joined_lineage,pedigree=INphylo$pedigree,nitt=13000,prior=prior1,verbose=F,burnin=3000)
            
            
            
            
            #Start filling in the data in the respective column
            drift_results[counter,1]<-tree_used
            drift_results[counter,2]<-repetition
            
            drift_results[counter,3]<-"positiveloss"
            
            if(length(unique(joined_ecology$subsistence))==2) drift_results[counter,4]<-summary(MCMC_regression_ecology)$solutions[2,5]
            
            if(length(unique(joined_alternative$alternative))==2) drift_results[counter,5]<-summary(MCMC_regression_alternative)$solutions[2,5]
            
            if(length(unique(joined_lineage$cladeA))==2) drift_results[counter,6]<-summary(MCMC_regression_lineage)$solutions[2,5]
            
            drift_results[counter,7]<-var(simulatedtips_neutral)
            
            drift_results[counter,8]<-median(joined_ecology[joined_ecology$subsistence==0,]$simulated)
            drift_results[counter,9]<-median(joined_ecology[joined_ecology$subsistence==1,]$simulated)
            
            drift_results[counter,10]<-median(joined_alternative[joined_alternative$alternative==0,]$simulated)
            drift_results[counter,11]<-median(joined_alternative[joined_alternative$alternative==1,]$simulated)
            
            drift_results[counter,12]<-median(joined_lineage[joined_lineage$cladeA==0,]$simulated)
            drift_results[counter,13]<-median(joined_lineage[joined_lineage$cladeA==1,]$simulated)
            
            drift_results[counter,14]<-predictor_variant_used
          
          counter<-counter+1
          
        
        
        
        
        
        
        
        
        # Loss according to ecoregion (50% of those living in forests or hunting/gathering)
        
        
        # Subset the tree and the data to assume that 50% of data is randomly missing
        
        
        
          
          randomdeletion<-function(x) (x*rbinom(1,1,0.5))
          Ecologytips<-sapply(as.numeric(Ecologymember),randomdeletion)
          
          
          if(Option=="WNAI") Ecologylosstips<-simulatedtips_neutral[Ecologytips==0]
          if(Option=="PamaNyungan") Ecologylosstips<-simulatedtips_neutral[Ecologytips==1]
          
          Ecologylosttree<-keep.tip(currenttree,names(Ecologylosstips))
          
          
          
          
          joined_ecology<-data.frame(as.numeric(Ecologylosstips),Ecologymember[names(Ecologymember) %in% names(Ecologylosstips)],names(Ecologylosstips))
          colnames(joined_ecology)<-c("simulated","subsistence","animal")
          
          
          INphylo<-inverseA(Ecologylosttree)
          MCMC_regression_ecology<-MCMCglmm(simulated ~ subsistence,random=~animal,data=joined_ecology,pedigree=INphylo$pedigree,nitt=13000,prior=prior1,verbose=F,burnin=3000)
          
          joined_alternative<-data.frame(as.numeric(Ecologylosstips),Alternativemember[names(Alternativemember) %in% names(Ecologylosstips)],names(Ecologylosstips))
          colnames(joined_alternative)<-c("simulated","alternative","animal")
          
          
          INphylo<-inverseA(Ecologylosttree)
          MCMC_regression_ecology<-MCMCglmm(simulated ~ alternative,random=~animal,data=joined_alternative,pedigree=INphylo$pedigree,nitt=13000,prior=prior1,verbose=F,burnin=3000)
          
          
          
          joined_lineage<-data.frame(as.numeric(Ecologylosstips),cladeAmember[names(cladeAmember) %in% names(Ecologylosstips)],names(Ecologylosstips))
          colnames(joined_lineage)<-c("simulated","cladeA","animal")
          
          
          MCMC_regression_lineage<-MCMCglmm(simulated ~ cladeA,random=~animal,data=joined_lineage,pedigree=INphylo$pedigree,nitt=13000,prior=prior1,verbose=F,burnin=3000)
          
          
          
          
          #Start filling in the data in the respective column
          drift_results[counter,1]<-tree_used
          drift_results[counter,2]<-repetition
          
          drift_results[counter,3]<-"ecologyloss"
          
          if(length(unique(joined_ecology$subsistence))==2) drift_results[counter,4]<-summary(MCMC_regression_ecology)$solutions[2,5]
          
          if(length(unique(joined_alternative$alternative))==2) drift_results[counter,5]<-summary(MCMC_regression_alternative)$solutions[2,5]
          
          if(length(unique(joined_lineage$cladeA))==2) drift_results[counter,6]<-summary(MCMC_regression_lineage)$solutions[2,5]
          
          drift_results[counter,7]<-var(simulatedtips_neutral)
          
          drift_results[counter,8]<-median(joined_ecology[joined_ecology$subsistence==0,]$simulated)
          drift_results[counter,9]<-median(joined_ecology[joined_ecology$subsistence==1,]$simulated)
          
          drift_results[counter,10]<-median(joined_alternative[joined_alternative$alternative==0,]$simulated)
          drift_results[counter,11]<-median(joined_alternative[joined_alternative$alternative==1,]$simulated)
          
          drift_results[counter,12]<-median(joined_lineage[joined_lineage$cladeA==0,]$simulated)
          drift_results[counter,13]<-median(joined_lineage[joined_lineage$cladeA==1,]$simulated)
          
          drift_results[counter,14]<-predictor_variant_used
          
          
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
        
          
          joined_ecology<-data.frame(as.numeric(horizontal10tips),Ecologymember,names(horizontal10tips))
          colnames(joined_ecology)<-c("simulated","subsistence","animal")
          
          
          INphylo<-inverseA(currenttree)
          MCMC_regression_ecology<-MCMCglmm(simulated ~ subsistence,random=~animal,data=joined_ecology,pedigree=INphylo$pedigree,nitt=13000,prior=prior1,verbose=F,burnin=3000)
          
          joined_alternative<-data.frame(as.numeric(horizontal10tips),Alternativemember,names(horizontal10tips))
          colnames(joined_alternative)<-c("simulated","alternative","animal")
          
          
          INphylo<-inverseA(currenttree)
          MCMC_regression_ecology<-MCMCglmm(simulated ~ alternative,random=~animal,data=joined_alternative,pedigree=INphylo$pedigree,nitt=13000,prior=prior1,verbose=F,burnin=3000)
          
          
          
          joined_lineage<-data.frame(as.numeric(horizontal10tips),cladeAmember,names(horizontal10tips))
          colnames(joined_lineage)<-c("simulated","cladeA","animal")
          
          
          MCMC_regression_lineage<-MCMCglmm(simulated ~ cladeA,random=~animal,data=joined_lineage,pedigree=INphylo$pedigree,nitt=13000,prior=prior1,verbose=F,burnin=3000)
          
          
          
          #Start filling in the data in the respective column
          drift_results[counter,1]<-tree_used
          drift_results[counter,2]<-repetition
          
          drift_results[counter,3]<-"horizontal10"
          
          if(length(unique(joined_ecology$subsistence))==2) drift_results[counter,4]<-summary(MCMC_regression_ecology)$solutions[2,5]
          
          if(length(unique(joined_alternative$alternative))==2) drift_results[counter,5]<-summary(MCMC_regression_alternative)$solutions[2,5]
          
          if(length(unique(joined_lineage$cladeA))==2) drift_results[counter,6]<-summary(MCMC_regression_lineage)$solutions[2,5]
          
          drift_results[counter,7]<-var(simulatedtips_neutral)
          
          drift_results[counter,8]<-median(joined_ecology[joined_ecology$subsistence==0,]$simulated)
          drift_results[counter,9]<-median(joined_ecology[joined_ecology$subsistence==1,]$simulated)
          
          drift_results[counter,10]<-median(joined_alternative[joined_alternative$alternative==0,]$simulated)
          drift_results[counter,11]<-median(joined_alternative[joined_alternative$alternative==1,]$simulated)
          
          drift_results[counter,12]<-median(joined_lineage[joined_lineage$cladeA==0,]$simulated)
          drift_results[counter,13]<-median(joined_lineage[joined_lineage$cladeA==1,]$simulated)
          
          drift_results[counter,14]<-predictor_variant_used
         
          
          counter<-counter+1
          

        
        
        
        
        
        # Horizontal transfer, 10% of societies change, adopt state from within their own clade
        
        
          
          societiestochange<-sample(simulatedtips_neutral[names(simulatedtips_neutral) %in% names(cladeAmember[cladeAmember==1])],round(0.1*nrow(data),0))
          societiestocopy<-sample(simulatedtips_neutral[names(simulatedtips_neutral) %in% names(cladeAmember[cladeAmember==1])],round(0.1*nrow(data),0))
          
          
          horizontal10tips<-simulatedtips_neutral
          
          for (i in 1:round(0.1*nrow(data),0)) {
            currentchange<-  societiestochange[i]
            horizontal10tips[names(horizontal10tips) %in% names(currentchange)]<-societiestocopy[i]
          }
          
          
          joined_ecology<-data.frame(as.numeric(horizontal10tips),Ecologymember,names(horizontal10tips))
          colnames(joined_ecology)<-c("simulated","subsistence","animal")
          
          
          INphylo<-inverseA(currenttree)
          MCMC_regression_ecology<-MCMCglmm(simulated ~ subsistence,random=~animal,data=joined_ecology,pedigree=INphylo$pedigree,nitt=13000,prior=prior1,verbose=F,burnin=3000)
          
          joined_alternative<-data.frame(as.numeric(horizontal10tips),Alternativemember,names(horizontal10tips))
          colnames(joined_alternative)<-c("simulated","alternative","animal")
          
          
          INphylo<-inverseA(currenttree)
          MCMC_regression_ecology<-MCMCglmm(simulated ~ alternative,random=~animal,data=joined_alternative,pedigree=INphylo$pedigree,nitt=13000,prior=prior1,verbose=F,burnin=3000)
          
          
          
          
          
          joined_lineage<-data.frame(as.numeric(horizontal10tips),cladeAmember,names(horizontal10tips))
          colnames(joined_lineage)<-c("simulated","cladeA","animal")
          
          
          MCMC_regression_lineage<-MCMCglmm(simulated ~ cladeA,random=~animal,data=joined_lineage,pedigree=INphylo$pedigree,nitt=13000,prior=prior1,verbose=F,burnin=3000)
          
          
          
          #Start filling in the data in the respective column
          drift_results[counter,1]<-tree_used
          drift_results[counter,2]<-repetition
          
          drift_results[counter,3]<-"horizontal10clade"
          
          if(length(unique(joined_ecology$subsistence))==2) drift_results[counter,4]<-summary(MCMC_regression_ecology)$solutions[2,5]
          
          if(length(unique(joined_alternative$alternative))==2) drift_results[counter,5]<-summary(MCMC_regression_alternative)$solutions[2,5]
          
          if(length(unique(joined_lineage$cladeA))==2) drift_results[counter,6]<-summary(MCMC_regression_lineage)$solutions[2,5]
          
          drift_results[counter,7]<-var(simulatedtips_neutral)
          
          drift_results[counter,8]<-median(joined_ecology[joined_ecology$subsistence==0,]$simulated)
          drift_results[counter,9]<-median(joined_ecology[joined_ecology$subsistence==1,]$simulated)
          
          drift_results[counter,10]<-median(joined_alternative[joined_alternative$alternative==0,]$simulated)
          drift_results[counter,11]<-median(joined_alternative[joined_alternative$alternative==1,]$simulated)
          
          drift_results[counter,12]<-median(joined_lineage[joined_lineage$cladeA==0,]$simulated)
          drift_results[counter,13]<-median(joined_lineage[joined_lineage$cladeA==1,]$simulated)
          
          drift_results[counter,14]<-predictor_variant_used
          
          counter<-counter+1
          

        
        
        # Horizontal transfer, 10% of societies change, adopt state from another society in same ecoregion
        
        
          
          societiestochange<-sample(simulatedtips_neutral[names(simulatedtips_neutral) %in% names(Ecologymember[Ecologymember==1])],round(0.1*nrow(data),0))
          societiestocopy<-sample(simulatedtips_neutral[names(simulatedtips_neutral) %in% names(Ecologymember[Ecologymember==1])],round(0.1*nrow(data),0))
          
          
          horizontal10tips<-simulatedtips_neutral
          
          for (i in 1:round(0.1*nrow(data),0)) {
            currentchange<-  societiestochange[i]
            horizontal10tips[names(horizontal10tips) %in% names(currentchange)]<-societiestocopy[i]
          }
          
          
          
          joined_ecology<-data.frame(as.numeric(horizontal10tips),Ecologymember,names(horizontal10tips))
          colnames(joined_ecology)<-c("simulated","subsistence","animal")
          
          
          INphylo<-inverseA(currenttree)
          MCMC_regression_ecology<-MCMCglmm(simulated ~ subsistence,random=~animal,data=joined_ecology,pedigree=INphylo$pedigree,nitt=13000,prior=prior1,verbose=F,burnin=3000)
          
          
          
          joined_alternative<-data.frame(as.numeric(horizontal10tips),Alternativemember,names(horizontal10tips))
          colnames(joined_alternative)<-c("simulated","alternative","animal")
          
          
          INphylo<-inverseA(currenttree)
          MCMC_regression_ecology<-MCMCglmm(simulated ~ alternative,random=~animal,data=joined_alternative,pedigree=INphylo$pedigree,nitt=13000,prior=prior1,verbose=F,burnin=3000)
          
          
          
          joined_lineage<-data.frame(as.numeric(horizontal10tips),cladeAmember,names(horizontal10tips))
          colnames(joined_lineage)<-c("simulated","cladeA","animal")
          
          
          MCMC_regression_lineage<-MCMCglmm(simulated ~ cladeA,random=~animal,data=joined_lineage,pedigree=INphylo$pedigree,nitt=13000,prior=prior1,verbose=F,burnin=3000)
          
          
          
          
          #Start filling in the data in the respective column
          drift_results[counter,1]<-tree_used
          drift_results[counter,2]<-repetition
          
          drift_results[counter,3]<-"horizontal10ecology"
          
          if(length(unique(joined_ecology$subsistence))==2) drift_results[counter,4]<-summary(MCMC_regression_ecology)$solutions[2,5]
          
          if(length(unique(joined_alternative$alternative))==2) drift_results[counter,5]<-summary(MCMC_regression_alternative)$solutions[2,5]
          
          if(length(unique(joined_lineage$cladeA))==2) drift_results[counter,6]<-summary(MCMC_regression_lineage)$solutions[2,5]
          
          drift_results[counter,7]<-var(simulatedtips_neutral)
          
          drift_results[counter,8]<-median(joined_ecology[joined_ecology$subsistence==0,]$simulated)
          drift_results[counter,9]<-median(joined_ecology[joined_ecology$subsistence==1,]$simulated)
          
          drift_results[counter,10]<-median(joined_alternative[joined_alternative$alternative==0,]$simulated)
          drift_results[counter,11]<-median(joined_alternative[joined_alternative$alternative==1,]$simulated)
          
          drift_results[counter,12]<-median(joined_lineage[joined_lineage$cladeA==0,]$simulated)
          drift_results[counter,13]<-median(joined_lineage[joined_lineage$cladeA==1,]$simulated)
          
          drift_results[counter,14]<-predictor_variant_used
          
          counter<-counter+1
          

          
          
          
          
          
          # Horizontal transfer, changes on 1% of branches / currently set to 5% - change the 20
          
          nodeswith_ht<-sample(names(internalnodes),round(nrow(data)/20,0))
          nodeswith_ht<-sort(nodeswith_ht)
          
          for (nodestochange in 1:round(nrow(data)/20,0)) {
          currentnode<-nodeswith_ht[nodestochange]
          currenttips <- tips(currenttree, currentnode)
          ht_tree<-keep.tip(currenttree,currenttips)  
  
          after_ht_tips<-fastBM(tree=ht_tree,a=sample(internalnodes[names(internalnodes) != currentnode],1),sig2=1)
         
          ht_tips<-simulatedtips_neutral
          ht_tips[names(ht_tips) %in% names(after_ht_tips)]<-after_ht_tips
          }
   

          
          joined_ecology<-data.frame(as.numeric(ht_tips),Ecologymember,names(ht_tips))
          colnames(joined_ecology)<-c("simulated","subsistence","animal")
          
          
          INphylo<-inverseA(currenttree)
          MCMC_regression_ecology<-MCMCglmm(simulated ~ subsistence,random=~animal,data=joined_ecology,pedigree=INphylo$pedigree,nitt=13000,prior=prior1,verbose=F,burnin=3000)
          
          
          joined_alternative<-data.frame(as.numeric(ht_tips),Alternativemember,names(ht_tips))
          colnames(joined_alternative)<-c("simulated","alternative","animal")
          
          
          INphylo<-inverseA(currenttree)
          MCMC_regression_ecology<-MCMCglmm(simulated ~ alternative,random=~animal,data=joined_alternative,pedigree=INphylo$pedigree,nitt=13000,prior=prior1,verbose=F,burnin=3000)
          
          
          
          joined_lineage<-data.frame(as.numeric(ht_tips),cladeAmember,names(ht_tips))
          colnames(joined_lineage)<-c("simulated","cladeA","animal")
          
          
          MCMC_regression_lineage<-MCMCglmm(simulated ~ cladeA,random=~animal,data=joined_lineage,pedigree=INphylo$pedigree,nitt=13000,prior=prior1,verbose=F,burnin=3000)
          
          
          
          
          #Start filling in the data in the respective column
          drift_results[counter,1]<-tree_used
          drift_results[counter,2]<-repetition
          
          drift_results[counter,3]<-"horizontalnodes"
          
          if(length(unique(joined_ecology$subsistence))==2) drift_results[counter,4]<-summary(MCMC_regression_ecology)$solutions[2,5]
          
          if(length(unique(joined_alternative$alternative))==2) drift_results[counter,5]<-summary(MCMC_regression_alternative)$solutions[2,5]
          
          if(length(unique(joined_lineage$cladeA))==2) drift_results[counter,6]<-summary(MCMC_regression_lineage)$solutions[2,5]
          
          drift_results[counter,7]<-var(simulatedtips_neutral)
          
          drift_results[counter,8]<-median(joined_ecology[joined_ecology$subsistence==0,]$simulated)
          drift_results[counter,9]<-median(joined_ecology[joined_ecology$subsistence==1,]$simulated)
          
          drift_results[counter,10]<-median(joined_alternative[joined_alternative$alternative==0,]$simulated)
          drift_results[counter,11]<-median(joined_alternative[joined_alternative$alternative==1,]$simulated)
          
          drift_results[counter,12]<-median(joined_lineage[joined_lineage$cladeA==0,]$simulated)
          drift_results[counter,13]<-median(joined_lineage[joined_lineage$cladeA==1,]$simulated)
          
          drift_results[counter,14]<-predictor_variant_used
          
          counter<-counter+1  
          
          
        
          # Horizontal transfer, changes on 1% of branches in dominant ecology / currently set to 5%
          
          Ecologytree<-keep.tip(currenttree,names(Ecologymember[Ecologymember ==1]))
          
          Ecologyreconstructed<-fastAnc(currenttree,Ecologymember)
          Ecologyreconstructed<-round(Ecologyreconstructed,0)
          
          
          nodeswith_ht<-sample(names(internalnodes[Ecologyreconstructed==1]),round(nrow(data)/20,0))
          nodeswith_ht<-sort(nodeswith_ht)
          
          
          for (nodestochange in 1:round(nrow(data)/20,0)) {
            currentnode<-nodeswith_ht[nodestochange]
            currenttips <- tips(currenttree, currentnode)
            ht_tree<-keep.tip(currenttree,currenttips)  
            
            after_ht_tips<-fastBM(tree=ht_tree,a=sample(internalnodes[Ecologyreconstructed==1],1),sig2=1)
            
            ht_tips<-simulatedtips_neutral
            ht_tips[names(ht_tips) %in% names(after_ht_tips)]<-after_ht_tips
          }
          
          
          
          joined_ecology<-data.frame(as.numeric(ht_tips),Ecologymember,names(ht_tips))
          colnames(joined_ecology)<-c("simulated","subsistence","animal")
          
          
          INphylo<-inverseA(currenttree)
          MCMC_regression_ecology<-MCMCglmm(simulated ~ subsistence,random=~animal,data=joined_ecology,pedigree=INphylo$pedigree,nitt=13000,prior=prior1,verbose=F,burnin=3000)
          
          
          joined_alternative<-data.frame(as.numeric(ht_tips),Alternativemember,names(ht_tips))
          colnames(joined_alternative)<-c("simulated","alternative","animal")
          
          
          INphylo<-inverseA(currenttree)
          MCMC_regression_ecology<-MCMCglmm(simulated ~ alternative,random=~animal,data=joined_alternative,pedigree=INphylo$pedigree,nitt=13000,prior=prior1,verbose=F,burnin=3000)
          
          
          
          joined_lineage<-data.frame(as.numeric(ht_tips),cladeAmember,names(ht_tips))
          colnames(joined_lineage)<-c("simulated","cladeA","animal")
          
          
          MCMC_regression_lineage<-MCMCglmm(simulated ~ cladeA,random=~animal,data=joined_lineage,pedigree=INphylo$pedigree,nitt=13000,prior=prior1,verbose=F,burnin=3000)
          
          
          
          
          #Start filling in the data in the respective column
          drift_results[counter,1]<-tree_used
          drift_results[counter,2]<-repetition
          
          drift_results[counter,3]<-"horizontalnodesecology"
          
          if(length(unique(joined_ecology$subsistence))==2) drift_results[counter,4]<-summary(MCMC_regression_ecology)$solutions[2,5]
          
          if(length(unique(joined_alternative$alternative))==2) drift_results[counter,5]<-summary(MCMC_regression_alternative)$solutions[2,5]
          
          if(length(unique(joined_lineage$cladeA))==2) drift_results[counter,6]<-summary(MCMC_regression_lineage)$solutions[2,5]
          
          drift_results[counter,7]<-var(simulatedtips_neutral)
          
          drift_results[counter,8]<-median(joined_ecology[joined_ecology$subsistence==0,]$simulated)
          drift_results[counter,9]<-median(joined_ecology[joined_ecology$subsistence==1,]$simulated)
          
          drift_results[counter,10]<-median(joined_alternative[joined_alternative$alternative==0,]$simulated)
          drift_results[counter,11]<-median(joined_alternative[joined_alternative$alternative==1,]$simulated)
          
          drift_results[counter,12]<-median(joined_lineage[joined_lineage$cladeA==0,]$simulated)
          drift_results[counter,13]<-median(joined_lineage[joined_lineage$cladeA==1,]$simulated)
          
          drift_results[counter,14]<-predictor_variant_used
          
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
          
          
          joined_ecology<-data.frame(as.numeric(horizontalNeighbor10tips),Ecologymember,names(horizontalNeighbor10tips))
          colnames(joined_ecology)<-c("simulated","subsistence","animal")
          
          
          INphylo<-inverseA(currenttree)
          MCMC_regression_ecology<-MCMCglmm(simulated ~ subsistence,random=~animal,data=joined_ecology,pedigree=INphylo$pedigree,nitt=13000,prior=prior1,verbose=F,burnin=3000)
          
          
          
          joined_alternative<-data.frame(as.numeric(horizontalNeighbor10tips),Alternativemember,names(horizontalNeighbor10tips))
          colnames(joined_alternative)<-c("simulated","alternative","animal")
          
          
          INphylo<-inverseA(currenttree)
          MCMC_regression_ecology<-MCMCglmm(simulated ~ alternative,random=~animal,data=joined_alternative,pedigree=INphylo$pedigree,nitt=13000,prior=prior1,verbose=F,burnin=3000)
          
          
          joined_lineage<-data.frame(as.numeric(horizontalNeighbor10tips),cladeAmember,names(horizontalNeighbor10tips))
          colnames(joined_lineage)<-c("simulated","cladeA","animal")
          
          
          MCMC_regression_lineage<-MCMCglmm(simulated ~ cladeA,random=~animal,data=joined_lineage,pedigree=INphylo$pedigree,nitt=13000,prior=prior1,verbose=F,burnin=3000)
          
          
          
          
          #Start filling in the data in the respective column
          drift_results[counter,1]<-tree_used
          drift_results[counter,2]<-repetition
          
          drift_results[counter,3]<-"horizontalneighbor"
          
          if(length(unique(joined_ecology$subsistence))==2) drift_results[counter,4]<-summary(MCMC_regression_ecology)$solutions[2,5]
          
          if(length(unique(joined_alternative$alternative))==2) drift_results[counter,5]<-summary(MCMC_regression_alternative)$solutions[2,5]
          
          if(length(unique(joined_lineage$cladeA))==2) drift_results[counter,6]<-summary(MCMC_regression_lineage)$solutions[2,5]
          
          drift_results[counter,7]<-var(simulatedtips_neutral)
          
          drift_results[counter,8]<-median(joined_ecology[joined_ecology$subsistence==0,]$simulated)
          drift_results[counter,9]<-median(joined_ecology[joined_ecology$subsistence==1,]$simulated)
          
          drift_results[counter,10]<-median(joined_alternative[joined_alternative$alternative==0,]$simulated)
          drift_results[counter,11]<-median(joined_alternative[joined_alternative$alternative==1,]$simulated)
          
          drift_results[counter,12]<-median(joined_lineage[joined_lineage$cladeA==0,]$simulated)
          drift_results[counter,13]<-median(joined_lineage[joined_lineage$cladeA==1,]$simulated)
          
          drift_results[counter,14]<-predictor_variant_used
          
          counter<-counter+1
          
          
          
        print("runningcounter")
        print(counter)
        print("repetition")
        print(repetition)
        print("currenttree")
        print(tree_used)
        
        
       
      
    } # end of the repetition loop
    
  if(Option=="WNAI") write.csv(drift_results[1:counter,],file="SimulateCulture_EcologyResults_Continuous_WNAI.csv")  
  if(Option=="PamaNyungan") write.csv(drift_results[1:counter,],file="SimulateCulture_EcologyResults_Continuous_PamaNyungan.csv")  
  
} #end of the tree variant loop
  
} #end of predictor variant loop
