## LDGM
Identifying Gene Regulatory Network Rewiring using Latent Differential Graphical Models.
![LDGM](/LDGM-workflow.png)


## A very brief introduction to LDGM

LDGM directly estimate the differential gene regulatory network between two biological conditions such as tumor subtypes. It has a clear advantage of utilizing much smaller sample size to achieve reliable differential network estimation compared with other Gaussian Graphical model based methods. To enable this, LDGM made several methodological contributions, such as   
- directly estimating the differential network from gene expression data without infering individual networks in each condition
- relaxing normality assumption on gene expression data

For more information including LDGM application on Breast cancer subtypes (Luminal A and Basel-like subtypes), please read our [paper](https://academic.oup.com/nar/article-abstract/44/17/e140/2468041) .

## Introduction to the repository

This repository contains:  

1. Matlab implementation of LDGM. The implementaion is written by Quanquan Gu (http://web.cs.ucla.edu/~qgu/). See folder "LDGM".  
2. LDGM application to generate a differential network. See folder "Stand_alone_example_by_LDGM". Running the example here to get started.

## Citation
To cite LDGM in your research work, please cite the corresponding paper  *Tian D, Gu Q and Ma J. 2016. Identifying gene regulatory network rewiring using latent differential graphical models. Nucleic Acids Research 44(17): e140.*

## Contact  
Dechao Tian at dechaotian@gmail.com.

Alternatively, please contact Jian Ma at jianma@cs.cmu.edu.

Computational Biology Department

School of Computer Science  

Carnegie Mellon University  

Please visit the [lab website](http://www.cs.cmu.edu/~jianma/) for more information about the Ma lab.

