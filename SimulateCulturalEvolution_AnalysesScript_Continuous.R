#------------------------------------------------------------------------------------------
# Analyses script of the simulation for the manuscript
#------------------------------------------------------------------------------------------

#Load necessary libraries

library(ape)
library(geiger)
library(phytools)
library(OUwie)
library(dplyr)
library(btw)

setwd("~/ownCloud/Documents/CulturalPhylogenetics/Code/PamaNyungan/")

WNAI_results<-read.csv("SimulateCulture_EcologyResults_Continuous_WNAI_1.csv")
WNAI_results_2<-read.csv("SimulateCulture_EcologyResults_Continuous_WNAI_v3.csv")

WNAI_allresults<-rbind(WNAI_results,WNAI_results_2)


WNAI_claderesults<-WNAI_allresults[WNAI_allresults$predictor_variable=="simmapattemptCladeA",]
WNAI_ecologyresults<-WNAI_allresults[WNAI_allresults$predictor_variable=="simmapattemptEcology",]
WNAI_alternativeresults<-WNAI_allresults[WNAI_allresults$predictor_variable=="simmapattemptAlternative",]


WNAI_notcladeresults<-WNAI_allresults[WNAI_allresults$predictor_variable!="simmapattemptCladeA",]
WNAI_notecologyresults<-WNAI_allresults[WNAI_allresults$predictor_variable!="simmapattemptEcology",]
WNAI_notalternativeresults<-WNAI_allresults[WNAI_allresults$predictor_variable!="simmapattemptAlternative",]



#WNAI false negative rate for clade
WNAI_claderesults %>% filter(pvalue_clade<0.05) %>% summarise(n_distinct(X))
WNAI_claderesults %>% filter(pvalue_clade>0.05) %>% summarise(n_distinct(X))
#0%

#WNAI false negative rate for forest
WNAI_ecologyresults %>% filter(pvalue_ecology<0.05) %>% summarise(n_distinct(X))
WNAI_ecologyresults %>% filter(pvalue_ecology>0.05) %>% summarise(n_distinct(X))
#1%

#WNAI false negative rate for alternative boats
WNAI_alternativeresults %>% filter(pvalue_alternative<0.05) %>% summarise(n_distinct(X))
WNAI_alternativeresults %>% filter(pvalue_alternative>0.05) %>% summarise(n_distinct(X))
#9%




#WNAI false positive rate for clade
WNAI_notcladeresults %>% filter(pvalue_clade<0.05) %>% summarise(n_distinct(X))
WNAI_notcladeresults %>% filter(pvalue_clade>0.05) %>% summarise(n_distinct(X))

(WNAI_notcladeresults %>% filter(pvalue_clade<0.05) %>% summarise(n_distinct(X))) / (  (WNAI_notcladeresults %>% filter(pvalue_clade>0.05) %>% summarise(n_distinct(X))) + (WNAI_notcladeresults %>% filter(pvalue_clade<0.05) %>% summarise(n_distinct(X))) )
# 6%


#WNAI false positive rate for ecology forest
WNAI_notecologyresults %>% filter(pvalue_ecology<0.05) %>% summarise(n_distinct(X))
WNAI_notecologyresults %>% filter(pvalue_ecology>0.05) %>% summarise(n_distinct(X))

(WNAI_notecologyresults %>% filter(pvalue_ecology<0.05) %>% summarise(n_distinct(X))) / (  (WNAI_notecologyresults %>% filter(pvalue_ecology>0.05) %>% summarise(n_distinct(X))) + (WNAI_notecologyresults %>% filter(pvalue_ecology<0.05) %>% summarise(n_distinct(X))) )
#15%

#WNAI false positive rate for alternative boats
WNAI_notalternativeresults %>% filter(pvalue_alternative<0.05) %>% summarise(n_distinct(X))
WNAI_notalternativeresults %>% filter(pvalue_alternative>0.05) %>% summarise(n_distinct(X))

(WNAI_notalternativeresults %>% filter(pvalue_alternative<0.05) %>% summarise(n_distinct(X))) / (  (WNAI_notalternativeresults %>% filter(pvalue_alternative>0.05) %>% summarise(n_distinct(X))) + (WNAI_notalternativeresults %>% filter(pvalue_alternative<0.05) %>% summarise(n_distinct(X))) )
# 8%




PamaNyungan_results<-read.csv("SimulateCulture_EcologyResults_Continuous_PamaNyungan_1.csv")
PamaNyungan_results_2<-read.csv("SimulateCulture_EcologyResults_Continuous_PamaNyungan_2.csv")
PamaNyungan_results_3<-read.csv("SimulateCulture_EcologyResults_Continuous_PamaNyungan_3.csv")
PamaNyungan_allresults<-PamaNyungan_results_3
PamaNyungan_allresults<-rbind(PamaNyungan_results_2,PamaNyungan_results_3)


PamaNyungan_claderesults<-PamaNyungan_allresults[PamaNyungan_allresults$predictor_variable=="simmapattemptCladeA",]
PamaNyungan_ecologyresults<-PamaNyungan_allresults[PamaNyungan_allresults$predictor_variable=="simmapattemptEcology",]
PamaNyungan_alternativeresults<-PamaNyungan_allresults[PamaNyungan_allresults$predictor_variable=="simmapattemptAlternative",]


PamaNyungan_notcladeresults<-PamaNyungan_allresults[PamaNyungan_allresults$predictor_variable!="simmapattemptCladeA",]
PamaNyungan_notecologyresults<-PamaNyungan_allresults[PamaNyungan_allresults$predictor_variable!="simmapattemptEcology",]
PamaNyungan_notalternativeresults<-PamaNyungan_allresults[PamaNyungan_allresults$predictor_variable!="simmapattemptAlternative",]



#PamaNyungan false negative rate for clade
PamaNyungan_claderesults %>% filter(pvalue_clade<0.05) %>% summarise(n_distinct(X))
PamaNyungan_claderesults %>% filter(pvalue_clade>0.05) %>% summarise(n_distinct(X))
#0%

#PamaNyungan false negative rate for ecology hunter gatherer
PamaNyungan_ecologyresults %>% filter(pvalue_clade<0.05) %>% summarise(n_distinct(X))
PamaNyungan_ecologyresults %>% filter(pvalue_clade>0.05) %>% summarise(n_distinct(X))
#28%

#PamaNyungan false negative rate for alternative red yellow
PamaNyungan_alternativeresults %>% filter(pvalue_alternative<0.05) %>% summarise(n_distinct(X))
PamaNyungan_alternativeresults %>% filter(pvalue_alternative>0.05) %>% summarise(n_distinct(X))
#45%




#PamaNyungan false positive rate for clade
PamaNyungan_notcladeresults %>% filter(pvalue_ecology<0.05) %>% summarise(n_distinct(X))
PamaNyungan_notcladeresults %>% filter(pvalue_ecology>0.05) %>% summarise(n_distinct(X))

(PamaNyungan_notcladeresults %>% filter(pvalue_ecology<0.05) %>% summarise(n_distinct(X))) / (  (PamaNyungan_notcladeresults %>% filter(pvalue_ecology>0.05) %>% summarise(n_distinct(X))) + (PamaNyungan_notcladeresults %>% filter(pvalue_ecology<0.05) %>% summarise(n_distinct(X))) )
# 33%


#PamaNyungan false positive rate for ecology hunter gatherer
PamaNyungan_notecologyresults %>% filter(pvalue_ecology<0.05) %>% summarise(n_distinct(X))
PamaNyungan_notecologyresults %>% filter(pvalue_ecology>0.05) %>% summarise(n_distinct(X))

(PamaNyungan_notecologyresults %>% filter(pvalue_ecology<0.05) %>% summarise(n_distinct(X))) / (  (PamaNyungan_notecologyresults %>% filter(pvalue_ecology>0.05) %>% summarise(n_distinct(X))) + (PamaNyungan_notecologyresults %>% filter(pvalue_ecology<0.05) %>% summarise(n_distinct(X))) )
#51%

#PamaNyungan false positive rate for alternative red/yellow
PamaNyungan_notalternativeresults %>% filter(pvalue_alternative<0.05) %>% summarise(n_distinct(X))
PamaNyungan_notalternativeresults %>% filter(pvalue_alternative>0.05) %>% summarise(n_distinct(X))

(PamaNyungan_notalternativeresults %>% filter(pvalue_alternative<0.05) %>% summarise(n_distinct(X))) / (  (PamaNyungan_notalternativeresults %>% filter(pvalue_alternative>0.05) %>% summarise(n_distinct(X))) + (PamaNyungan_notalternativeresults %>% filter(pvalue_alternative<0.05) %>% summarise(n_distinct(X))) )
# 42%








#PamaNyungan false positive rate for ecology hunter gatherer in non-random loss
(PamaNyungan_notecologyresults %>% filter(pvalue_ecology>0.05,Sample=="positiveloss") %>% summarise(n_distinct(X))) / (  (PamaNyungan_notecologyresults %>% filter(pvalue_ecology>0.05,Sample=="positiveloss") %>% summarise(n_distinct(X))) + (PamaNyungan_notecologyresults %>% filter(pvalue_ecology<0.05,Sample=="positiveloss") %>% summarise(n_distinct(X))) )
#PamaNyungan false negative rate for ecology hunter gatherer in non-random loss
(PamaNyungan_ecologyresults %>% filter(pvalue_ecology>0.05,Sample=="positiveloss") %>% summarise(n_distinct(X))) / (  (PamaNyungan_ecologyresults %>% filter(pvalue_ecology>0.05,Sample=="positiveloss") %>% summarise(n_distinct(X))) + (PamaNyungan_ecologyresults %>% filter(pvalue_ecology<0.05,Sample=="positiveloss") %>% summarise(n_distinct(X))) )


(WNAI_notecologyresults %>% filter(pvalue_ecology<0.05,Sample=="positiveloss") %>% summarise(n_distinct(X))) / (  (WNAI_notecologyresults %>% filter(pvalue_ecology>0.05,Sample=="positiveloss") %>% summarise(n_distinct(X))) + (WNAI_notecologyresults %>% filter(pvalue_ecology<0.05,Sample=="positiveloss") %>% summarise(n_distinct(X))) )
(WNAI_ecologyresults %>% filter(pvalue_ecology>0.05,Sample=="positiveloss") %>% summarise(n_distinct(X))) / (  (WNAI_ecologyresults %>% filter(pvalue_ecology>0.05,Sample=="positiveloss") %>% summarise(n_distinct(X))) + (WNAI_ecologyresults %>% filter(pvalue_ecology<0.05,Sample=="positiveloss") %>% summarise(n_distinct(X))) )

#PamaNyungan false positive rate for clade in non-random loss
(PamaNyungan_notcladeresults %>% filter(pvalue_clade>0.05,Sample=="positiveloss") %>% summarise(n_distinct(X))) / (  (PamaNyungan_notcladeresults %>% filter(pvalue_clade>0.05,Sample=="positiveloss") %>% summarise(n_distinct(X))) + (PamaNyungan_notcladeresults %>% filter(pvalue_clade<0.05,Sample=="positiveloss") %>% summarise(n_distinct(X))) )
#PamaNyungan false negative rate for clade in non-random loss
(PamaNyungan_claderesults %>% filter(pvalue_clade<0.05,Sample=="positiveloss") %>% summarise(n_distinct(X))) / (  (PamaNyungan_claderesults %>% filter(pvalue_clade>0.05,Sample=="positiveloss") %>% summarise(n_distinct(X))) + (PamaNyungan_claderesults %>% filter(pvalue_clade<0.05,Sample=="positiveloss") %>% summarise(n_distinct(X))) )

(WNAI_notcladeresults %>% filter(pvalue_clade<0.05,Sample=="positiveloss") %>% summarise(n_distinct(X))) / (  (WNAI_notcladeresults %>% filter(pvalue_clade>0.05,Sample=="positiveloss") %>% summarise(n_distinct(X))) + (WNAI_notcladeresults %>% filter(pvalue_clade<0.05,Sample=="positiveloss") %>% summarise(n_distinct(X))) )
(WNAI_claderesults %>% filter(pvalue_clade>0.05,Sample=="positiveloss") %>% summarise(n_distinct(X))) / (  (WNAI_claderesults %>% filter(pvalue_clade>0.05,Sample=="positiveloss") %>% summarise(n_distinct(X))) + (WNAI_claderesults %>% filter(pvalue_clade<0.05,Sample=="positiveloss") %>% summarise(n_distinct(X))) )



#PamaNyungan false positive rate for ecology hunter gatherer in non-random loss
(PamaNyungan_notecologyresults %>% filter(pvalue_ecology>0.05,Sample=="horizontal10") %>% summarise(n_distinct(X))) / (  (PamaNyungan_notecologyresults %>% filter(pvalue_ecology>0.05,Sample=="horizontal10") %>% summarise(n_distinct(X))) + (PamaNyungan_notecologyresults %>% filter(pvalue_ecology<0.05,Sample=="horizontal10") %>% summarise(n_distinct(X))) )
#PamaNyungan false negative rate for ecology hunter gatherer in non-random loss
(PamaNyungan_ecologyresults %>% filter(pvalue_ecology>0.05,Sample=="horizontal10") %>% summarise(n_distinct(X))) / (  (PamaNyungan_ecologyresults %>% filter(pvalue_ecology>0.05,Sample=="horizontal10") %>% summarise(n_distinct(X))) + (PamaNyungan_ecologyresults %>% filter(pvalue_ecology<0.05,Sample=="horizontal10") %>% summarise(n_distinct(X))) )


(WNAI_notecologyresults %>% filter(pvalue_ecology<0.05,Sample=="horizontal10") %>% summarise(n_distinct(X))) / (  (WNAI_notecologyresults %>% filter(pvalue_ecology>0.05,Sample=="horizontal10") %>% summarise(n_distinct(X))) + (WNAI_notecologyresults %>% filter(pvalue_ecology<0.05,Sample=="horizontal10") %>% summarise(n_distinct(X))) )
(WNAI_ecologyresults %>% filter(pvalue_ecology>0.05,Sample=="horizontal10") %>% summarise(n_distinct(X))) / (  (WNAI_ecologyresults %>% filter(pvalue_ecology>0.05,Sample=="horizontal10") %>% summarise(n_distinct(X))) + (WNAI_ecologyresults %>% filter(pvalue_ecology<0.05,Sample=="horizontal10") %>% summarise(n_distinct(X))) )

#PamaNyungan false positive rate for clade in non-random loss
(PamaNyungan_notcladeresults %>% filter(pvalue_clade>0.05,Sample=="horizontal10") %>% summarise(n_distinct(X))) / (  (PamaNyungan_notcladeresults %>% filter(pvalue_clade>0.05,Sample=="horizontal10") %>% summarise(n_distinct(X))) + (PamaNyungan_notcladeresults %>% filter(pvalue_clade<0.05,Sample=="horizontal10") %>% summarise(n_distinct(X))) )
#PamaNyungan false negative rate for clade in non-random loss
(PamaNyungan_claderesults %>% filter(pvalue_clade<0.05,Sample=="horizontal10") %>% summarise(n_distinct(X))) / (  (PamaNyungan_claderesults %>% filter(pvalue_clade>0.05,Sample=="horizontal10") %>% summarise(n_distinct(X))) + (PamaNyungan_claderesults %>% filter(pvalue_clade<0.05,Sample=="horizontal10") %>% summarise(n_distinct(X))) )

(WNAI_notcladeresults %>% filter(pvalue_clade<0.05,Sample=="horizontal10") %>% summarise(n_distinct(X))) / (  (WNAI_notcladeresults %>% filter(pvalue_clade>0.05,Sample=="horizontal10") %>% summarise(n_distinct(X))) + (WNAI_notcladeresults %>% filter(pvalue_clade<0.05,Sample=="horizontal10") %>% summarise(n_distinct(X))) )
(WNAI_claderesults %>% filter(pvalue_clade>0.05,Sample=="horizontal10") %>% summarise(n_distinct(X))) / (  (WNAI_claderesults %>% filter(pvalue_clade>0.05,Sample=="horizontal10") %>% summarise(n_distinct(X))) + (WNAI_claderesults %>% filter(pvalue_clade<0.05,Sample=="horizontal10") %>% summarise(n_distinct(X))) )






#PamaNyungan false positive rate for ecology hunter gatherer in non-random loss
(PamaNyungan_notecologyresults %>% filter(pvalue_ecology>0.05,Tree=="Latetree") %>% summarise(n_distinct(X))) / (  (PamaNyungan_notecologyresults %>% filter(pvalue_ecology>0.05,Tree=="Latetree") %>% summarise(n_distinct(X))) + (PamaNyungan_notecologyresults %>% filter(pvalue_ecology<0.05,Tree=="Latetree") %>% summarise(n_distinct(X))) )
#PamaNyungan false negative rate for ecology hunter gatherer in non-random loss
(PamaNyungan_ecologyresults %>% filter(pvalue_ecology>0.05,Tree=="Latetree") %>% summarise(n_distinct(X))) / (  (PamaNyungan_ecologyresults %>% filter(pvalue_ecology>0.05,Tree=="Latetree") %>% summarise(n_distinct(X))) + (PamaNyungan_ecologyresults %>% filter(pvalue_ecology<0.05,Tree=="Latetree") %>% summarise(n_distinct(X))) )


(WNAI_notecologyresults %>% filter(pvalue_ecology<0.05,Tree=="Latetree") %>% summarise(n_distinct(X))) / (  (WNAI_notecologyresults %>% filter(pvalue_ecology>0.05,Tree=="Latetree") %>% summarise(n_distinct(X))) + (WNAI_notecologyresults %>% filter(pvalue_ecology<0.05,Tree=="Latetree") %>% summarise(n_distinct(X))) )
(WNAI_ecologyresults %>% filter(pvalue_ecology>0.05,Tree=="Latetree") %>% summarise(n_distinct(X))) / (  (WNAI_ecologyresults %>% filter(pvalue_ecology>0.05,Tree=="Latetree") %>% summarise(n_distinct(X))) + (WNAI_ecologyresults %>% filter(pvalue_ecology<0.05,Tree=="Latetree") %>% summarise(n_distinct(X))) )

#PamaNyungan false positive rate for clade in non-random loss
(PamaNyungan_notcladeresults %>% filter(pvalue_clade<0.05,Tree=="Latetree") %>% summarise(n_distinct(X))) / (  (PamaNyungan_notcladeresults %>% filter(pvalue_clade>0.05,Tree=="Latetree") %>% summarise(n_distinct(X))) + (PamaNyungan_notcladeresults %>% filter(pvalue_clade<0.05,Tree=="Latetree") %>% summarise(n_distinct(X))) )
#PamaNyungan false negative rate for clade in non-random loss
(PamaNyungan_claderesults %>% filter(pvalue_clade>0.05,Tree=="Latetree") %>% summarise(n_distinct(X))) / (  (PamaNyungan_claderesults %>% filter(pvalue_clade>0.05,Tree=="Latetree") %>% summarise(n_distinct(X))) + (PamaNyungan_claderesults %>% filter(pvalue_clade<0.05,Tree=="Latetree") %>% summarise(n_distinct(X))) )

(WNAI_notcladeresults %>% filter(pvalue_clade<0.05,Tree=="Latetree") %>% summarise(n_distinct(X))) / (  (WNAI_notcladeresults %>% filter(pvalue_clade>0.05,Tree=="Latetree") %>% summarise(n_distinct(X))) + (WNAI_notcladeresults %>% filter(pvalue_clade<0.05,Tree=="Latetree") %>% summarise(n_distinct(X))) )
(WNAI_claderesults %>% filter(pvalue_clade>0.05,Tree=="Latetree") %>% summarise(n_distinct(X))) / (  (WNAI_claderesults %>% filter(pvalue_clade>0.05,Tree=="Latetree") %>% summarise(n_distinct(X))) + (WNAI_claderesults %>% filter(pvalue_clade<0.05,Tree=="Latetree") %>% summarise(n_distinct(X))) )









WNAI_directional<-read.csv("SimulateCulture_DirectionalResults_Continuous_WNAI.csv")
WNAI_nodirectional<-read.csv("SimulateCulture_NoDirectionalResults_Continuous_WNAI.csv")

PamaNy_directional<-read.csv("SimulateCulture_DirectionalResults_Continuous_PamaNyungan.csv")
PamaNy_nodirectional<-read.csv("SimulateCulture_NoDirectionalResults_Continuous_PamaNyungan.csv")



correct_WNAI_directional_sample<-(WNAI_directional %>% filter(pvalue_trend<0.05) %>% group_by(Sample) %>% summarise(n_distinct(X))) 
false_WNAI_directional_sample<-(WNAI_directional %>% filter(pvalue_trend>0.05) %>% group_by(Sample) %>% summarise(n_distinct(X))) 

correct_WNAI_nodirectional_sample<-(WNAI_nodirectional %>% filter(pvalue_trend>0.05) %>% group_by(Sample) %>% summarise(n_distinct(X))) 
false_WNAI_nodirectional_sample<-(WNAI_nodirectional %>% filter(pvalue_trend<0.05) %>% group_by(Sample) %>% summarise(n_distinct(X))) 


correct_PamaNy_directional_sample<-(PamaNy_directional %>% filter(pvalue_trend<0.05) %>% group_by(Sample) %>% summarise(n_distinct(X))) 
false_PamaNy_directional_sample<-(PamaNy_directional %>% filter(pvalue_trend>0.05) %>% group_by(Sample) %>% summarise(n_distinct(X))) 

correct_PamaNy_nodirectional_sample<-(PamaNy_nodirectional %>% filter(pvalue_trend>0.05) %>% group_by(Sample) %>% summarise(n_distinct(X))) 
false_PamaNy_nodirectional_sample<-(PamaNy_nodirectional %>% filter(pvalue_trend<0.05) %>% group_by(Sample) %>% summarise(n_distinct(X))) 





correct_WNAI_directional_tree<-(WNAI_directional %>% filter(pvalue_trend<0.05) %>% group_by(Tree) %>% summarise(n_distinct(X))) 
false_WNAI_directional_tree<-(WNAI_directional %>% filter(pvalue_trend>0.05) %>% group_by(Tree) %>% summarise(n_distinct(X))) 

correct_WNAI_nodirectional_tree<-(WNAI_nodirectional %>% filter(pvalue_trend>0.05) %>% group_by(Tree) %>% summarise(n_distinct(X))) 
false_WNAI_nodirectional_tree<-(WNAI_nodirectional %>% filter(pvalue_trend<0.05) %>% group_by(Tree) %>% summarise(n_distinct(X))) 


correct_PamaNy_directional_tree<-(PamaNy_directional %>% filter(pvalue_trend<0.05) %>% group_by(Tree) %>% summarise(n_distinct(X))) 
false_PamaNy_directional_tree<-(PamaNy_directional %>% filter(pvalue_trend>0.05) %>% group_by(Tree) %>% summarise(n_distinct(X))) 

correct_PamaNy_nodirectional_tree<-(PamaNy_nodirectional %>% filter(pvalue_trend>0.05) %>% group_by(Tree) %>% summarise(n_distinct(X))) 
false_PamaNy_nodirectional_tree<-(PamaNy_nodirectional %>% filter(pvalue_trend<0.05) %>% group_by(Tree) %>% summarise(n_distinct(X))) 






drift_results %>% filter(chisq_p_selectionwronglysupported<0.05,horizontaltransmission %in% "No") %>% summarise(n_distinct(X))

drift_results %>% filter(chisq_p_selectionwronglysupported<0.05,horizontaltransmission %in% "No") %>% group_by(RateOfChange) %>% summarise(n_distinct(X))
drift_results %>% filter(horizontaltransmission %in% "No") %>% group_by(RateOfChange) %>% summarise(n_distinct(X))


drift_results %>% filter(chisq_p_selectionwronglysupported<0.05,horizontaltransmission %in% "No") %>% group_by(NumberOfVariantsObserved) %>% summarise(n_distinct(X))
drift_results %>% filter(horizontaltransmission %in% "No") %>% group_by(NumberOfVariantsObserved) %>% summarise(n_distinct(X))


drift_results %>% filter(chisq_p_selectionwronglysupported<0.05,horizontaltransmission %in% "Yes") %>% group_by(NumberOfVariantsObserved) %>% summarise(n_distinct(X))
drift_results %>% filter(horizontaltransmission %in% "Yes") %>% group_by(NumberOfVariantsObserved) %>% summarise(n_distinct(X))

drift_results %>% filter(chisq_p_selectionwronglysupported<0.05,horizontaltransmission %in% "Yes") %>% group_by(RateOfChange) %>% summarise(n_distinct(X))
drift_results %>% filter(horizontaltransmission %in% "Yes") %>% group_by(RateOfChange) %>% summarise(n_distinct(X))


drift_results %>% filter(chisq_p_selectionwronglysupported<0.05,horizontaltransmission %in% "No") %>% group_by(Tree) %>% summarise(n_distinct(X))
drift_results %>% filter(horizontaltransmission %in% "No") %>% group_by(Tree) %>% summarise(n_distinct(X))


drift_results %>% filter(chisq_p_selectionwronglysupported<0.05,horizontaltransmission %in% "Yes") %>% group_by(Tree) %>% summarise(n_distinct(X))
drift_results %>% filter(horizontaltransmission %in% "Yes") %>% group_by(Tree) %>% summarise(n_distinct(X))



nonconverged<-drift_results %>%  filter(BT_Multistate_nonconvergence %in% "TRUE") %>% group_by(Sample) %>% summarise(n_distinct(X))
converged<-drift_results %>%  filter(BT_Multistate_nonconvergence %in% "FALSE") %>% group_by(Sample) %>% summarise(n_distinct(X))


drift_results %>% filter(chisq_BT_MultitstateFree_better>0.05,horizontaltransmission %in% "Yes") %>% group_by(Sample) %>% summarise(n_distinct(X))
drift_results %>% filter(chisq_BT_MultitstateFree_better<0.05,horizontaltransmission %in% "Yes") %>% group_by(Sample) %>% summarise(n_distinct(X))

drift_results %>% filter(chisq_BT_MultitstateFree_better>0.05,horizontaltransmission %in% "No") %>% group_by(NumberOfVariantsObserved) %>% summarise(n_distinct(X))
drift_results %>% filter(chisq_BT_MultitstateFree_better<0.05,horizontaltransmission %in% "No") %>% group_by(NumberOfVariantsObserved) %>% summarise(n_distinct(X))


drift_results %>% filter(chisq_BT_MultitstateFree_better>0.05,horizontaltransmission %in% "Yes",BT_Multistate_nonconvergence %in% "FALSE") %>% group_by(Sample) %>% summarise(n_distinct(X))
drift_results %>% filter(chisq_BT_MultitstateFree_better<0.05,horizontaltransmission %in% "Yes",BT_Multistate_nonconvergence %in% "FALSE") %>% group_by(Sample) %>% summarise(n_distinct(X))


drift_results %>% filter(chisq_BT_MultitstateFree_better>0.05,horizontaltransmission %in% "No",BT_Multistate_nonconvergence %in% "FALSE") %>% group_by(NumberOfVariantsObserved) %>% summarise(n_distinct(X))
drift_results %>% filter(chisq_BT_MultitstateFree_better<0.05,horizontaltransmission %in% "No",BT_Multistate_nonconvergence %in% "FALSE") %>% group_by(NumberOfVariantsObserved) %>% summarise(n_distinct(X))

drift_results %>% filter(chisq_BT_MultitstateFree_better>0.05, BT_Multistate_nonconvergence %in% "FALSE") %>% group_by(Tree) %>% summarise(n_distinct(X))
drift_results %>% filter(chisq_BT_MultitstateFree_better<0.05, BT_Multistate_nonconvergence %in% "FALSE") %>% group_by(Tree) %>% summarise(n_distinct(X))


drift_results %>% filter(chisq_BT_CladeAdependent_better>0.05) %>% group_by(Tree) %>% summarise(n_distinct(X))
drift_results %>% filter(chisq_BT_CladeAdependent_better<0.05) %>% group_by(Tree) %>% summarise(n_distinct(X))

drift_results %>% filter(chisq_BT_CladeAdependent_better<0.05) %>% group_by(chisq_BT_CladeAdependent_better<0.05) %>% summarise(mean(phylosigK))


drift_results %>% filter(chisq_BT_Forestdependent_better>0.05) %>% group_by(Tree) %>% summarise(n_distinct(X))
drift_results %>% filter(chisq_BT_Forestdependent_better<0.05) %>% group_by(Tree) %>% summarise(n_distinct(X))


drift_results %>% filter(chisq_BT_Forestdependent_better<0.05) %>% group_by(Sample) %>% summarise(n_distinct(X))
drift_results %>% filter(chisq_BT_Forestdependent_better>0.05) %>% group_by(Sample) %>% summarise(n_distinct(X))


drift_results %>%  group_by(chisq_BT_MultitstateFree_better<0.5) %>% summarise(mean(phylosigK))




selection_results %>% filter(horizontaltransmission %in% "No") %>% summarise(n_distinct(X))

selection_results %>% filter(chisq_p_selectioncorrectlysupported>0.05,horizontaltransmission %in% "No") %>% summarise(n_distinct(X))

selection_results %>% filter(chisq_p_selectioncorrectlysupported>0.05,horizontaltransmission %in% "No") %>% group_by(RateOfChange) %>% summarise(n_distinct(X))
selection_results %>% filter(horizontaltransmission %in% "No") %>% group_by(RateOfChange) %>% summarise(n_distinct(X))


selection_results %>% filter(chisq_p_selectioncorrectlysupported>0.05,horizontaltransmission %in% "No") %>% group_by(NumberOfVariantsObserved) %>% summarise(n_distinct(X))
selection_results %>% filter(horizontaltransmission %in% "No") %>% group_by(NumberOfVariantsObserved) %>% summarise(n_distinct(X))


selection_results %>% filter(chisq_p_selectioncorrectlysupported>0.05,horizontaltransmission %in% "Yes") %>% group_by(NumberOfVariantsObserved) %>% summarise(n_distinct(X))
selection_results %>% filter(horizontaltransmission %in% "Yes") %>% group_by(NumberOfVariantsObserved) %>% summarise(n_distinct(X))

selection_results %>% filter(chisq_p_selectioncorrectlysupported>0.05,horizontaltransmission %in% "Yes") %>% group_by(RateOfChange) %>% summarise(n_distinct(X))
selection_results %>% filter(horizontaltransmission %in% "Yes") %>% group_by(RateOfChange) %>% summarise(n_distinct(X))


selection_results %>% filter(chisq_p_selectioncorrectlysupported>0.05,horizontaltransmission %in% "No") %>% group_by(Tree) %>% summarise(n_distinct(X))
selection_results %>% filter(horizontaltransmission %in% "No") %>% group_by(Tree) %>% summarise(n_distinct(X))


selection_results %>% filter(chisq_p_selectioncorrectlysupported>0.05,horizontaltransmission %in% "Yes") %>% group_by(Tree) %>% summarise(n_distinct(X))
selection_results %>% filter(horizontaltransmission %in% "Yes") %>% group_by(Tree) %>% summarise(n_distinct(X))


selection_results %>% filter(chisq_BT_CladeAdependent_better<0.05,horizontaltransmission %in% "No") %>% group_by(Tree) %>% summarise(n_distinct(X))
selection_results %>% filter(chisq_BT_CladeAdependent_better>0.05,horizontaltransmission %in% "No") %>% group_by(Tree) %>% summarise(n_distinct(X))


selection_results %>% filter(chisq_BT_CladeAdependent_better<0.05,horizontaltransmission %in% "Yes") %>% group_by(Tree) %>% summarise(n_distinct(X))
selection_results %>% filter(chisq_BT_CladeAdependent_better>0.05,horizontaltransmission %in% "Yes") %>% group_by(Tree) %>% summarise(n_distinct(X))


selection_results %>% filter(chisq_BT_Forestdependent_better<0.05,horizontaltransmission %in% "No") %>% group_by(Tree) %>% summarise(n_distinct(X))
selection_results %>% filter(chisq_BT_Forestdependent_better>0.05,horizontaltransmission %in% "No") %>% group_by(Tree) %>% summarise(n_distinct(X))


drift_results %>% filter(chisq_BT_CladeAdependent_better<0.05,horizontaltransmission %in% "Yes") %>% group_by(Tree) %>% summarise(n_distinct(X))
drift_results %>% filter(chisq_BT_CladeAdependent_better>0.05,horizontaltransmission %in% "Yes") %>% group_by(Tree) %>% summarise(n_distinct(X))


drift_results %>% filter(chisq_BT_Forestdependent_better<0.05,horizontaltransmission %in% "No") %>% group_by(Tree) %>% summarise(n_distinct(X))
drift_results %>% filter(chisq_BT_Forestdependent_better>0.05,horizontaltransmission %in% "No") %>% group_by(Tree) %>% summarise(n_distinct(X))



selection_results %>% filter(chisq_BT_CladeAdependent_better<0.05) %>% group_by(Sample) %>% summarise(n_distinct(X))
selection_results %>% filter(chisq_BT_CladeAdependent_better>0.05) %>% group_by(Sample) %>% summarise(n_distinct(X))

drift_results %>% filter(chisq_BT_CladeAdependent_better<0.05) %>% group_by(Sample) %>% summarise(n_distinct(X))
drift_results %>% filter(chisq_BT_CladeAdependent_better>0.05) %>% group_by(Sample) %>% summarise(n_distinct(X))





plot(samplesize_analysis[1:12,]$proportion_drift_selectionwronglysupported~samplesize_analysis[1:12,]$samplesize,xlim=c(0,200),ylim=c(0,1))
points(samplesize_analysis$proportion_selection_selectionwronglynotsupported~samplesize_analysis$samplesize,col="red")
points(samplesize_analysis$drift_diversityrepresented~samplesize_analysis$samplesize,col="black",pch=8)
points(samplesize_analysis$selection_diversityrepresented~samplesize_analysis$samplesize,col="red",pch=8)




#outcome variables: 
#drift: proportion wrongly assigned as selection (~20%) - false positives
#drift: proportion of total variants recovered
#drift: Clade A Bayestraits wrongly supported
#drift: Forest Bayestraits wrongly supported
#drift: Bayestraits Multistate dependent wrongly supported

#selection: proportion wrongly not detected as selection (~50%) - false negatives
#selection: proportion of total variants recovered
#selection: Clade A Bayestraits wrongly supported
#selection: Forest Bayestraits wrongly supported
#selection: Bayestraits Multistate dependent wrongly supported

#predictor variables:
#sample size
#rate of change
#number of variants in model
#horizontal transmission
#tree shape
#for selection: to one or pathway



drift_results_filtered<-drift_results[ complete.cases(drift_results$chisq_p_selectionwronglysupported) , ]


dat_list <- list(
  chisq_p_selectionwronglysupported = drift_results_filtered$chisq_p_selectionwronglysupported,
  samplesize = standardize(drift_results_filtered$samplesize),
  horizontaltransmission = as.integer(drift_results_filtered$horizontaltransmission)
)

m8.1 <- quap( alist(
  chisq_p_selectionwronglysupported ~ dnorm( mu , sigma ) , 
  mu <- a + b*samplesize , 
  a ~ dnorm( 1 , 1 ) ,
  b ~ dnorm( 0 , 1 ) ,
  sigma ~ dexp( 1 ) ) , data=dat_list )
precis(m8.1)


set.seed(7)
prior <- extract.prior( m8.1 )
# set up the plot dimensions
plot( NULL , xlim=c(0,1) , ylim=c(0.5,1.5) ,
      xlab="ruggedness" , ylab="log GDP" ) abline( h=min(dd$log_gdp_std) , lty=2 ) abline( h=max(dd$log_gdp_std) , lty=2 )
# draw 50 lines from the prior
rugged_seq <- seq( from=-0.1 , to=1.1 , length.out=30 )
mu <- link( m8.1 , post=prior , data=data.frame(rugged_std=rugged_seq) ) for ( i in 1:50 ) lines( rugged_seq , mu[i,] , col=col.alpha("black",0.3) )


Latitudedistance<-as.matrix(dist(WNAIlocations$Latitude))
Longitudedistance<-as.matrix(dist(WNAIlocations$Longitude))
Totaldistance<-Latitudedistance+Longitudedistance
diag(Totaldistance)<-0
colnames(Totaldistance)<-row.names(WNAIlocations)
row.names(Totaldistance)<-row.names(WNAIlocations)
distancetree = nj(Totaldistance)
plot.phylo(distancetree,edge.width=3,label.offset=1, cex=0.4)



plot.phylo(Grafentree,type="fan",edge.color = (simmapattemptEcology$mapped.edge[,1]>0.5)+1,show.tip.label = F,edge.width=3)

plot.phylo(Onetree,edge.width=3,label.offset=1, cex=0.4,direction="leftwards")

