############################################
### Run commands listed bleow sequentially
############################################

time python Scripts/BRCA/BRCA.py

matlab -singleCompThread -nojvm -nodesktop -nodisplay < Scripts/BRCA//DGM/DGM_comparison.m
