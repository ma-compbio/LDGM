#########################################
#### Here are codes for 
#### reproducing validation results
#### when p = 30 and m = 2.
#### Run the codes sequentially  
##########################################

############################
#### Generate data, including
#### individual networks,
#### Differential network \Delta
#### Correlations matrices
############################
1) Differential_netwk_generation_rho.py

############################
#### Estimate differetil 
#### networks by DGM
############################
2) DGM/DGM_comparison.m

############################
#### Estimate differetil 
#### networks by Glasso
############################
3) Glasso_comparison.m

############################
#### Compute AUC under ROC curves
#### Draw boxplots of AUC
############################
4) ROC_slim.py