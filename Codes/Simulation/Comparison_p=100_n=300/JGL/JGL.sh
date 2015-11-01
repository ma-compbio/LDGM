#!/bin/sh 
cd "$(dirname "$0")" 
Rscript JGL_fused_10.R 
Rscript JGL_fused_1.R 
Rscript JGL_fused_0.1.R 
Rscript JGL_fused_0.01.R 
Rscript JGL_fused_0.001.R 
Rscript JGL_fused_0.0001.R 