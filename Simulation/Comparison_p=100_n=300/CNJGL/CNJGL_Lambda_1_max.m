function [lambda_1_max] = Find_Lambda_max(S1, S2, lambda_1, lambda_2, n1, n2)
	[Theta_1,Theta_2] = ADM_CNJGL(S1,S2,lambda_1,lambda_2, n1, n2);
	% Theta_1 = abs(Theta_1) > 10e-5;
	% Theta_2 = abs(Theta_2) > 10e-5;
	Theta_1( abs(Theta_1) < 10e-5 ) = 0;
	Theta_2( abs(Theta_2) < 10e-5 ) = 0;
	hdelta = Theta_1 - Theta_2;

	% ind = ~ any(hdelta(:));
	iter = 0;
	if nnz(triu(hdelta, 1)) == 0
		while nnz(triu(hdelta, 1)) == 0 & iter < 10
			lambda_1 = lambda_1 / 10;
			[Theta_1,Theta_2] = ADM_CNJGL(S1,S2,lambda_1,lambda_2, n1, n2);
			Theta_1( abs(Theta_1) < 10e-5 ) = 0;
			Theta_2( abs(Theta_2) < 10e-5 ) = 0;
			% Theta_1 = abs(Theta_1) > 10e-5;
			% Theta_2 = abs(Theta_2) > 10e-5;
			hdelta = Theta_1 - Theta_2;
			iter = iter +1;
		end
		lambda_1_max = lambda_1 * 10;
	end

	if nnz(triu(hdelta, 1)) ~= 0
		while nnz(triu(hdelta, 1)) ~= 0 & iter < 10
			lambda_1 = lambda_1 * 10;
			[Theta_1,Theta_2] = ADM_CNJGL(S1,S2,lambda_1,lambda_2, n1, n2);
			Theta_1( abs(Theta_1) < 10e-5 ) = 0;
			Theta_2( abs(Theta_2) < 10e-5 ) = 0;
			% Theta_1 = abs(Theta_1) > 10e-5;
			% Theta_2 = abs(Theta_2) > 10e-5;
			hdelta = Theta_1 - Theta_2;
			iter = iter +1;
		end
		lambda_1_max = lambda_1;
	end
end
