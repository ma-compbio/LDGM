function [lambda_1_min] = CNJGL_Lambda_1_min(S1, S2, lambda_1_max, n_target, ... 
															lambda_2, n1, n2)
	lambda_1 = lambda_1_max;
	
	tp = 0;
	iter = 0;
	while tp < n_target & iter < 10
		lambda_1 = 0.5 * lambda_1;
		[Theta_1,Theta_2] = ADM_CNJGL(S1,S2,lambda_1, lambda_2, n1, n2);
		% Theta_1 = abs(Theta_1) > 10e-5;
		% Theta_2 = abs(Theta_2) > 10e-5;
		Theta_1( abs(Theta_1) < 10e-5) = 0;
		Theta_2( abs(Theta_2) < 10e-5) = 0;
		hdelta = Theta_1 - Theta_2;
		tp = nnz(triu(hdelta, 1));
		iter = iter + 1;
	end
	lambda_1_min = lambda_1;
end

