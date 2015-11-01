clc
clear
close all
format long g

cd('/home/dctian/GGM')

% ############################
p = 30;
n_simu = 500;
num_diff_edge = 2;
Tissue_key = {'Lung', 'Whole Blood', 'Brain', 'Heart'};
% ###### 0-base
i = 2;
j = 1;
% ##############################
k = 30;


DGGM_name = 'DGGM';
folder = strcat(Tissue_key{i + 1}, '_vs_', Tissue_key{j+1}, '_p=', num2str(p), '_nd=', num2str(num_diff_edge));
	
addpath(strcat('Scripts/Real_Data_Validation/', folder,'/', DGGM_name))

folder = strcat('Results/Real_Data_Validation/', folder, '/');

file_output = strcat(folder, DGGM_name, '_delta_vs_hdelta.txt');
if exist(file_output, 'file') == 2
	delete(file_output)
end 
file_output2 = strcat(folder, DGGM_name, '_delta_union_vs_hdelta.txt');
if exist(file_output2, 'file') == 2
	delete(file_output2)
end 

lambda_output = strcat(folder, DGGM_name, '_lambdas.txt');
if exist(lambda_output, 'file') == 2
	delete(lambda_output)
end	

iter = 0;
parameter_combin = strcat(num2str(p), ...
				'_nd=', num2str(num_diff_edge), ... 
				'_iter=',num2str(iter));
file_hsigma1 = strcat(folder, 'hsigma_1_p=', parameter_combin,'.txt')
file_hsigma2 = strcat(folder,'hsigma_2_p=', parameter_combin,'.txt');
file_delta = strcat(folder,'delta_p=', parameter_combin,'.txt');
file_delta_union = strcat(folder,'delta_union_p=', parameter_combin,'.txt');

hsigma1 = csvread(file_hsigma1);
hsigma2 = csvread(file_hsigma2);

delta = csvread(file_delta);
delta = delta + 1;
num_de = size(delta, 1);

delta_union = csvread(file_delta_union);
delta_union = delta_union + 1;
num_de_union = size(delta_union, 1);
n_target = p * (p-1) / 2;

lambda_max = Lambda_max(hsigma1, hsigma2);
lambda_min = Lambda_min(hsigma1, hsigma2, lambda_max, n_target);
lambda_seq = linspace(lambda_min, lambda_max, k+1);
lambda_seq = fliplr(lambda_seq(1:end-1));
dlmwrite(lambda_output, lambda_seq, '-append');
Result = Estimate_delta(hsigma1, hsigma2, lambda_seq, delta, delta_union);
result = Result(:, 1:2);
result = result(:);
result(end + 1, 1) = num_de;
result = result.';
dlmwrite(file_output, result, '-append')

result_union = Result(:, 3:4);
result_union = result_union(:);
result_union(end + 1, 1) = num_de_union;
result_union = result_union.';
dlmwrite(file_output2, result_union, '-append')
for iter = 1:(n_simu-1)
	parameter_combin = strcat(num2str(p), ...
					'_nd=', num2str(num_diff_edge), ... 
					'_iter=',num2str(iter));
	file_hsigma1 = strcat(folder, 'hsigma_1_p=', parameter_combin,'.txt');
	file_hsigma2 = strcat(folder,'hsigma_2_p=', parameter_combin,'.txt');
	file_delta = strcat(folder,'delta_p=', parameter_combin,'.txt');
	file_delta_union = strcat(folder,'delta_union_p=', parameter_combin,'.txt');

	hsigma1 = csvread(file_hsigma1);
	hsigma2 = csvread(file_hsigma2);

	delta = csvread(file_delta);
	delta = delta + 1;
	num_de = size(delta, 1);

	delta_union = csvread(file_delta_union);
	delta_union = delta_union + 1;
	num_de_union = size(delta_union, 1);

	Result = Estimate_delta(hsigma1, hsigma2, lambda_seq, delta, delta_union);
	result = Result(:, 1:2);
	result = result(:);
	result(end + 1, 1) = num_de;
	result = result.';
	dlmwrite(file_output, result, '-append')

	result_union = Result(:, 3:4);
	result_union = result_union(:);
	result_union(end + 1, 1) = num_de_union;
	result_union = result_union.';
	dlmwrite(file_output2, result_union, '-append')
end	
exit;
