clc
clear
close all
format long g

cd('/home/dctian/GGM')

k = 30;


DGGM_name = 'DGGM';
folder = 'BRCA';
	
addpath(strcat('Scripts/', folder,'/', DGGM_name))

folder = strcat('Results/', folder, '/');

lambda_output = strcat(folder, DGGM_name, '_lambdas.txt');
if exist(lambda_output, 'file') == 2
	delete(lambda_output)
end	

file_hsigma1 = strcat(folder, 'hsigma_1.txt');
file_hsigma2 = strcat(folder,'hsigma_2.txt');

hsigma1 = csvread(file_hsigma1);
hsigma2 = csvread(file_hsigma2);

p = size(hsigma1, 1);
n_target = p * (p-1) / 2;

lambda_max = Lambda_max(hsigma1, hsigma2);
lambda_min = Lambda_min(hsigma1, hsigma2, lambda_max, n_target);
lambda_seq = linspace(lambda_min, lambda_max, k+1);
lambda_seq = fliplr(lambda_seq(1:end-1))
dlmwrite(lambda_output, lambda_seq, '-append');

for i = 1:30
	hdelta = differential_graph(hsigma1, hsigma2, lambda_seq(i));
	hdelta = triu(hdelta, 1);
	if nnz(hdelta) > 0
		[row, col, v] = find(hdelta);
		edgelist = [row, col];
		% 1-base to 0-base
		edgelist = edgelist - 1;
		file_edge = strcat(folder,'Edgelist_i=', num2str(i), '.txt');
		dlmwrite(file_edge, edgelist,'\t')
	end	
end
exit;
