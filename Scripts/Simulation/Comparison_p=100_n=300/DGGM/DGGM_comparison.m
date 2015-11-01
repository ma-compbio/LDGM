clc
clear
close all
format long g

cd('/home/dctian/GGM')

DGGM_name = 'DGGM';

folder = 'Comparison_p=100_n=300';

addpath(strcat('Scripts/Simulation/', folder,'/', DGGM_name))


p = 100;
rho_all = [0.05, 0.1, 0.2];
n = 300;
k = 30;
n_simu = 30;


folder = strcat('Results/Simulation/', folder, '/');

input_file_m = strcat(folder, 'Parameters/m.txt');
m_all = csvread(input_file_m);

file_output = strcat(folder, DGGM_name, '_delta_vs_hdelta.txt');
if exist(file_output, 'file') == 2
	delete(file_output)
end 

lambda_output = strcat(folder, DGGM_name, '_lambdas.txt');
if exist(lambda_output, 'file') == 2
	delete(lambda_output)
end	

for m=m_all
	for rho=rho_all
		iter = 0;
		parameter_combin = strcat(num2str(p), ...
						'_m=', num2str(m), '_rho1=', num2str(rho), ... 
						'_iter=',num2str(iter));
		file_hsigma1 = strcat(folder, 'hsigma_1_p=', parameter_combin,'.txt');
		file_hsigma2 = strcat(folder,'hsigma_2_p=', parameter_combin,'.txt');
		file_delta = strcat(folder,'delta_p=', parameter_combin,'.txt');
		file_ROC = strcat(folder, DGGM_name,'_p=', parameter_combin);

		hsigma1 = csvread(file_hsigma1);
		hsigma2 = csvread(file_hsigma2);

		delta = csvread(file_delta);
		num_de = nnz(triu(delta, 1));
		n_target = p * (p-1) / 2;

		lambda_max = Lambda_max(hsigma1, hsigma2);
		lambda_min = Lambda_min(hsigma1, hsigma2, lambda_max, n_target);
		lambda_seq = linspace(lambda_min, lambda_max, k+1);
		lambda_seq = fliplr(lambda_seq(1:end-1));
		dlmwrite(lambda_output, lambda_seq, '-append');

		result = Estimate_delta(hsigma1, hsigma2, lambda_seq, delta, file_ROC);
		result = result(:);
		result(end + 1, 1) = num_de;
		result = result.';
		dlmwrite(file_output, result, '-append')
		for iter = 1:(n_simu-1)
			parameter_combin = strcat(num2str(p), ...
							'_m=', num2str(m), '_rho1=', num2str(rho), ... 
							'_iter=',num2str(iter));
			file_hsigma1 = strcat(folder, 'hsigma_1_p=', parameter_combin,'.txt');
			file_hsigma2 = strcat(folder,'hsigma_2_p=', parameter_combin,'.txt');
			file_delta = strcat(folder,'delta_p=', parameter_combin,'.txt');
			file_ROC = strcat(folder, DGGM_name,'_p=', parameter_combin);

			hsigma1 = csvread(file_hsigma1);
			hsigma2 = csvread(file_hsigma2);
			delta = csvread(file_delta);
			num_de = nnz(triu(delta, 1));
			
			result = Estimate_delta(hsigma1, hsigma2, lambda_seq, delta, file_ROC);
			result = result(:);
			result(end + 1, 1) = num_de;
			result = result.';
			dlmwrite(file_output, result, '-append')
		end	
	end	
end
exit;

