clc
clear
close all
format long g

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% First change the working directory
%%% for example by 
%%% cd('/Stand_alone_example_by_DGM')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath('LDGM')

file_hsigma1 = 'Data/hsigma_1.txt';
file_hsigma2 = 'Data/hsigma_2.txt';
file_TFs = 'Data/TFs.txt'

hsigma1 = csvread(file_hsigma1);
hsigma2 = csvread(file_hsigma2);

lambda = 0.26924
hdelta = differential_graph(hsigma1, hsigma2, lambda);
hdelta = triu(hdelta, 1);
[row, col, v] = find(hdelta);
edgelist = [row, col];

%%%%% generate edgelist where nodes are gene symbol
fid = fopen(file_TFs, 'r');
formatSpec = '%s';
N = 4;
C_text = textscan(fid,formatSpec,N,'Delimiter','\t');
C_data0 = textscan(fid,'%d %s %s %s');
fclose(fid)
% 0-base to 1-base
Gene_id = C_data0{1} + 1;
Gene_name = C_data0{2};
Gene_name = cellstr(Gene_name);

edgelist_gene = num2cell(edgelist);
mask = ismember(edgelist, Gene_id);
edgelist_gene(mask) = Gene_name(edgelist(mask))
