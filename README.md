# CulturalMacroevolution_Simulation

This repository contains the code and data related to the manuscript Lukas, Towner & Borgerhoff-Mulder "The Potential to Infer the Historical Pattern of Cultural Macroevolution as illustrated by the Western North American Indian Database". 

The purpose of these files is to provide a setup to simulate the evolution of cultural traits across phylogenies to understand the potential risks of false inferences. 

In the manuscript, we present results from two example datasets: the WNAI societies and Pama Nyungan societies. For these, the input data are in the "example data" folder (for each dataset giving the phylogeny, geographic locations, and ecological categorisation of the the societies). 

To run the simulations in R, use the "SimulateCulturalEvolution_MainScript.R" file. The main script calls other code files from the folder "code". As explained in the main script, you might want to refer to these other files directly in case you want to repeat the simulations across your own dataset.

The particular output of the simulations we analysed for the manuscript are provided as tables in the folder "Results". These can be loaded into R and comparisons of which models fit best for a given simulation can, for example, be assessed using a likelihood ratio test (e.g. for a simulation with four variant 1-pchisq(2*(log likelihood of a reconstruction with a model assuming that some transitions are more likely than others minus log likelihood of a reconstruction with a model assuming that all transitions are equally likely), 11 ) where 11 reflects the degree of freedoms, based on the difference between the 12 different transition rates between all four variants in the free model and the 1 transition rate in the equal rates model. 
