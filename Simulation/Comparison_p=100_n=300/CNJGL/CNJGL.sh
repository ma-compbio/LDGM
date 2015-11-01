#!/bin/sh 
cd "$(dirname "$0")" 
module load matlab 
#matlab -singleCompThread -nodesktop -nodisplay < CNJGL_lambda_2=10n.m 
matlab -singleCompThread -nodesktop -nodisplay < CNJGL_lambda_2=1n.m 
matlab -singleCompThread -nodesktop -nodisplay < CNJGL_lambda_2=0.1n.m 
matlab -singleCompThread -nodesktop -nodisplay < CNJGL_lambda_2=0.01n.m 
matlab -singleCompThread -nodesktop -nodisplay < CNJGL_lambda_2=0.001n.m 
matlab -singleCompThread -nodesktop -nodisplay < CNJGL_lambda_2=0.0001n.m 
