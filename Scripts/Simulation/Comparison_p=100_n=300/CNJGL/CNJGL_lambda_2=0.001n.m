clc
clear
close all
format long g

cd('/home/dctian/GGM')

CNJGL_name = 'CNJGL';

folder = 'Comparison_p=100_n=300';

addpath(strcat('Scripts/Simulation/', folder,'/', CNJGL_name))

sender = 'stark198507@gmail.com';
psswd = 'stark1985';

setpref('Internet','E_mail',sender);
setpref('Internet','SMTP_Server','smtp.gmail.com');
setpref('Internet','SMTP_Username',sender);
setpref('Internet','SMTP_Password',psswd);

props = java.lang.System.getProperties;
props.setProperty('mail.smtp.auth','true');
props.setProperty('mail.smtp.socketFactory.class', ...
                  'javax.net.ssl.SSLSocketFactory');
props.setProperty('mail.smtp.socketFactory.port','465');

recipient = 'dechaotian@gmail.com';

p = 100;
rho_all = [0.05, 0.1, 0.2];
n = 300;
k = 30;
n_simu = 30;

lambda_1 = 0.1 * n;
lambda_2_sample = 0.001;
lambda_2 = lambda_2_sample * n;

n1 = n;
n2 = n;


folder = strcat('Results/Simulation/', folder, '/');

input_file_m = strcat(folder, 'Parameters/m.txt');
m_all = csvread(input_file_m);

file_output = strcat(folder, CNJGL_name, '_lambda2=',num2str(lambda_2_sample),'n_delta_vs_hdelta.txt');
if exist(file_output, 'file') == 2
	delete(file_output)
end 

lambda_output = strcat(folder, CNJGL_name, '_lambda2=',num2str(lambda_2_sample), 'n_lambdas.txt');
if exist(lambda_output, 'file') == 2
	delete(lambda_output)
end	


for m=m_all
	for rho=rho_all
		
		iter = 0;
		file_hsigma1 = strcat(folder, 'hsigma_1_p=', num2str(p), ...
						'_m=', num2str(m), '_rho1=', num2str(rho), ... 
						'_iter=',num2str(iter),'.txt');
		file_hsigma2 = strcat(folder,'hsigma_2_p=', num2str(p), ... 
						'_m=', num2str(m), '_rho1=', num2str(rho), ... 
						'_iter=', num2str(iter),'.txt');
		file_delta = strcat(folder,'delta_p=', num2str(p),'_m=', num2str(m), ...
						'_rho1=',num2str(rho),'_iter=',num2str(iter),'.txt');
		file_ROC = strcat(folder, CNJGL_name, '_p=', num2str(p),'_m=', num2str(m), ...
						'_rho1=',num2str(rho),'_iter=',num2str(iter), ... 
						'_lambda2=', num2str(lambda_2_sample), 'n');
						
		hsigma1 = csvread(file_hsigma1);
		hsigma2 = csvread(file_hsigma2);
		delta = csvread(file_delta);

		num_de = nnz(triu(delta, 1));
		n_target = p * (p - 1) / 2;

		lambda_1_max = CNJGL_Lambda_1_max(hsigma1, hsigma2, lambda_1, ... 
											lambda_2, n1, n2);
		lambda_1_min = CNJGL_Lambda_1_min(hsigma1, hsigma2, ...
							lambda_1_max, n_target, lambda_2, n1, n2);
		lambda_1_seq = linspace(lambda_1_min, lambda_1_max, k + 1);
		lambda_1_seq = fliplr(lambda_1_seq(1:end-1));
		dlmwrite(lambda_output, lambda_1_seq, '-append');

		result = CNJGL_estimate_delta(hsigma1, hsigma2, lambda_1_seq, ...
											lambda_2, n1, n2, delta, file_ROC);
		result = result(:);
		result(end + 1, 1) = num_de;
		result = result.';
		dlmwrite(file_output, result, '-append')
		for iter = 1:(n_simu-1)
			file_hsigma1 = strcat(folder, 'hsigma_1_p=', num2str(p), ...
							'_m=', num2str(m), '_rho1=', num2str(rho), ... 
							'_iter=',num2str(iter),'.txt');
			file_hsigma2 = strcat(folder,'hsigma_2_p=', num2str(p), ... 
							'_m=', num2str(m), '_rho1=', num2str(rho), ... 
							'_iter=', num2str(iter),'.txt');
			file_delta = strcat(folder,'delta_p=', num2str(p),'_m=', num2str(m), ...
							'_rho1=',num2str(rho),'_iter=',num2str(iter),'.txt');
			file_ROC = strcat(folder, CNJGL_name, '_p=', num2str(p),'_m=', num2str(m), ...
							'_rho1=',num2str(rho),'_iter=',num2str(iter), ...
							'_lambda2=', num2str(lambda_2_sample), 'n');

			hsigma1 = csvread(file_hsigma1);
			hsigma2 = csvread(file_hsigma2);
			delta = csvread(file_delta);
			num_de = nnz(triu(delta, 1));
			
			result = CNJGL_estimate_delta(hsigma1, hsigma2, lambda_1_seq, ...
												lambda_2, n1, n2, delta, file_ROC);
			result = result(:);
			result(end + 1, 1) = num_de;
			result = result.';
			dlmwrite(file_output, result, '-append')
		end	
	end
end
subject = sprintf('CNJGL_lambda2=%dn with p=%d is done', lambda_2_sample, p );
message = 'Congratulations';
sendmail(recipient, subject, message);


exit;

