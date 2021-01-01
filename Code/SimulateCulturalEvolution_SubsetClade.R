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
