# QBI 2025 
## PETase Variant Km Predictors 


PETase is an enzyme found in Ideonella sakaiensis bacteria. Its main function is to digest polyethylene terephthalate (PET), plastic commonly used to package food and beverages. 

By optimizing PETase efficiency, we can increase its potential for industrial usage to degrade plastic and, therefore, reduce plastic pollution. 

Here, we trained XGBoost and Linear Regression models using data similar to enzyme sequence and Km data from UniProt. The model aims to predict Km values of mutated PETase sequences to reduce costs associated with the experimental derivation of Km values. 

We wrote a script to generate random mutants of the PETase amino acid sequence. Then, we ran the mutated sequences to predict their KM values. 

This repo contains the script to generate the mutated sequences XGBoost, and linear regression notebooks. 



