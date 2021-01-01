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
}
