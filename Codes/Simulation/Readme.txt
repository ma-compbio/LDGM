#########################################
#### Here are codes for 
#### reproducing simulation results
#### when p = 100 and n = 300.
#### Run the codes in the folder 
#### Comparison_p=100_n=300 sequentially  
##########################################

############################
#### Generate data, including
#### Scale-free networks,
#### Differential network \Delta
#### Correlations matrices
############################
1) SF_network.py 

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
#### Estimate differetil 
#### networks by JGL with
#### fixed \lambda_2
############################
4) JGL/JGL.sh

############################
#### Estimate differetil 
#### networks by CNJGL with
#### fixed \lambda_2
############################
5) CNJGL/CNJGL.sh

############################
#### Drop ROC curve, 
#### Precision-Recall curve
#### Compute AUC
############################
6) ROC_PR.py