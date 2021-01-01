# These simulations and their outcomes are described in:
# Lukas, D., Towner, M., & Mulder, M. B. (2020). The Potential to Infer the Historical Pattern of Cultural Macroevolution as Illustrated by the Western North American Indian Societies.
# https://osf.io/preprints/socarxiv/tjvgy/

# Source file to load data for the examples

# Load relevant phylogenies from GitHub


if(Option=="WNAI") {
  
  #Load WNAI tree constructed from the language classifications
  Americantree<-read.nexus(url("https://raw.githubusercontent.com/dieterlukas/CulturalMacroevolution_Simulation/master/ExampleData/WNAI_tree_forsimulation.nex"))
}



if(Option=="PamaNyungan") {
  
  #Load Pama Nyungan tree provided by Bouckaert et al. 2018 https://doi.org/10.1038/s41559-018-0489-3
  PamaNyungantree<-read.nexus(url("https://raw.githubusercontent.com/dieterlukas/CulturalMacroevolution_Simulation/master/ExampleData/PamaNyungan_tree_forsimulation.nex"))
}


#------------------------------------------------------------------------------------------
# Load associated data from GitHub

# For both sets of societies, there is information on one ecological variable. 

# For the WNAI, the variable classifies the ecoregion that each society predominantly uses (32 different ecosystems; the data are from a version of the WNAI dataset provided by Dow & Eff http://intersci.ss.uci.edu/wiki/index.php/Materials_for_cross-cultural_research)
if(Option=="WNAI") {
  WNAIdata<-read_csv(url("https://raw.githubusercontent.com/dieterlukas/CulturalMacroevolution_Simulation/master/ExampleData/WNAI_data_forsimulation.csv"))
  WNAIdata<-data.frame(WNAIdata)
  
  # In addition, we load the geographic location for each society (latitude/longitude) for the horizontal transmission among neighbors; the data are from the same dataset
  WNAIlocations<-read_csv(url("https://raw.githubusercontent.com/dieterlukas/CulturalMacroevolution_Simulation/master/ExampleData/WNAI_locations_forsimulation.csv"))
  WNAIlocations<-data.frame(WNAIlocations)
  
}

# For the Pama Nyungan, the variable classifies the main mode of subsistence (hunter-gatherer versus food produce; the data are from Derungs et al. 2018 https://github.com/curdon/linguisticDensity_ProcB_derungsEtAl)
if(Option=="PamaNyungan") {
  PamaNyungandata<-read_csv(url("https://raw.githubusercontent.com/dieterlukas/CulturalMacroevolution_Simulation/master/ExampleData/PamaNyungan_data_forsimulation.csv"))
  PamaNyungandata<-data.frame(PamaNyungandata)
  rownames(PamaNyungandata)<-PamaNyungandata$society
  
  # In addition, we load the geographic location for each society (latitude/longitude) for the horizontal transmission among neighbors; the data are from the publication that describes the phylogeny, based on the Bowern, Claire (2016). Chirila: Contemporary and Historical Resources for the Indigenous Languages of Australia. Language Documentation and Conservation. Vol 10. http://www.pamanyungan.net/chirila/
  PamaNyunganlocations<-read_csv(url("https://raw.githubusercontent.com/dieterlukas/CulturalMacroevolution_Simulation/master/ExampleData/PamaNyungan_locations_forsimulation.csv"))
  PamaNyunganlocations<-data.frame(PamaNyunganlocations)
  rownames(PamaNyunganlocations)<-PamaNyunganlocations$Language
}

