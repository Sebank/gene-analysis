# gene-analysis
Code for project and master thesis

analysis pipeline.R

Showcases one way to utilize limma-voom and DESeq2 to perform gene alalysis. 
A simple GLMM has also been implemented. 
Does preprocessinig of data, removing data thought to be disruptive for the rest. 
Shows results and comparisons between them and does analysis of the different fits. 

gene wise simulation.R

Simulation of negative binomial distribution with relevant size factor estimated from count data in total.R. 
Compares estimated distribution of GLMM and GLM for different correlation and proportion of correlation when the data is generated as GLMM. 
