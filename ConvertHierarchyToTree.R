# The first part defines functions that will be used for constructing the trees

## recursion function
traverse <- function(a,i,innerl){
  if(i < (ncol(df))){
    alevelinner <- as.character(unique(df[which(as.character(df[,i])==a),i+1]))
    desc <- NULL
    if(length(alevelinner) == 1) (newickout <- traverse(alevelinner,i+1,innerl))
    else {
      for(b in alevelinner) desc <- c(desc,traverse(b,i+1,innerl))
      il <- NULL; if(innerl==TRUE) il <- a
      (newickout <- paste("(",paste(desc,collapse=","),")",il,sep=""))
    }
  }
  else { (newickout <- a) }
}

## data.frame to newick function
df2newick <- function(df, innerlabel=FALSE){
  alevel <- as.character(unique(df[,1]))
  newick <- NULL
  for(x in alevel) newick <- c(newick,traverse(x,1,innerlabel))
  (newick <- paste("(",paste(newick,collapse=","),");",sep=""))
}

df <- data.frame(x=c('A','A','B','B','B'), y=c('Ab','Ac','Ba', 'Ba','Bd'), z=c('Abb','Acc','Bad', 'Bae','Bdd'))
myNewick <- df2newick(df)


# Load necessary libraries to maninpulate data and the phylogeny
library(ape)
library(tidyr)
library(dplyr)


# Read the file with the hierarchical classification of societies across nine levels of language family affiliation
americantribes<-read.csv("https://raw.githubusercontent.com/dieterlukas/CulturalMacroevolution_Simulation/master/WNAI_languageclassification_numerical.csv")
americantribes<-americantribes[,2:10]
df<-americantribes
AmericanNewick <- df2newick(americantribes)


Americantree <- read.tree(text=AmericanNewick)
plot(Americantree)

# Read the file with the original classification and the original labels for each society
americantribes<-read.csv("https://raw.githubusercontent.com/dieterlukas/CulturalMacroevolution_Simulation/master/WNAI_languageclassification_labels.csv")
americantribes<-americantribes[,2:10]
for (i in 1:ncol(americantribes)) {
  americantribes[,i]<-as.numeric(americantribes[,i])
  }

americantribesjoined<-americantribes
taxonomiclevels<-colnames(americantribes)
for (i in 2:ncol(americantribes)) {
  americantribesjoined[,i]<-paste(americantribesjoined[,i-1],americantribesjoined[,i],sep="")
}

df<-americantribesjoined
AmericanNewick <- df2newick(americantribesjoined)
Americantree <- read.tree(text=AmericanNewick)
Americantreelabels<-matrix(Americantree$tip.label,ncol=1)
reclassify<-matrix(c(americantribes$L9,americantribesjoined$L9),nrow=172,ncol=2)
colnames(reclassify)<-c("labels","L9")
colnames(Americantreelabels)<-"L9"
Americantreelabels<-left_join(as.data.frame(Americantreelabels), as.data.frame(reclassify), by = "L9")
Americantree$tip.label<-as.character(Americantreelabels$labels)
plot(Americantree)
write.nexus(Americantree,file="americantree.nex")


