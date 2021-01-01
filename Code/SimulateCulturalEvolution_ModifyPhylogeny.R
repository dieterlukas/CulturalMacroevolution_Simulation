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
  
  PamaNyungantree<-force.ultrametric(PamaNyungantree)  
  
  #Add branch lengths to the tree based on Grafen's method (all tips equidistant from root, branch length depends onnumber of nodes between root and tip)
  Grafentree<-compute.brlen(PamaNyungantree,method="Grafen")
  
  #Add branch lengths to the tree assuming that all branches have the same length of one
  Onetree<-compute.brlen(PamaNyungantree,1)
  
  #Add branch lengths to the tree with early diversification and long branches to the tips
  Earlytree<-compute.brlen(PamaNyungantree,method="Grafen",power=0.25)
  
  #Add branch lengths to the tree with lots of recent diversification and long branches between clades
  Latetree<-compute.brlen(PamaNyungantree,method="Grafen",power=1.5)
  
  #Some analyses need a rooted, fully bifurcating tree
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
