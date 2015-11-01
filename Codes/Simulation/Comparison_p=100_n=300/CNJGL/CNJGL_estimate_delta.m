function result = CNJGL_estimate_delta(hsigma1, hsigma2, lambda_1_seq, lambda_2, n1, n2, delta, file_ROC)
	result = [];
	delta = delta(triu(true(size(delta)),1));
	for lambda_1 = lambda_1_seq
		% file_ROC_lambda = strcat(file_ROC, '_lambda=', num2str(lambda_1),'.txt');
		[Theta_1, Theta_2] = ADM_CNJGL(hsigma1, hsigma2,lambda_1,lambda_2,n1,n2);
		Theta_1( abs(Theta_1) < 10e-5) = 0;
		Theta_2( abs(Theta_2) < 10e-5) = 0;
		% Theta_1 = abs(Theta_1) > 10e-5;
		% Theta_2 = abs(Theta_2) > 10e-5;
		hdelta = Theta_1 - Theta_2;
		hdelta = hdelta(triu(true(size(hdelta)),1));

		% ROC_delta_hdelta(hdelta, delta, file_ROC_lambda);
		[x, y] = Compare_delta_hdelta(hdelta, delta);
		result = [result; [x, y]];
	end
end

function [x, y] = Compare_delta_hdelta(hdelta, delta)
	x = nnz(hdelta ~= 0);
	y = nnz(hdelta ~= 0 & delta ~= 0);
end

function ROC_delta_hdelta(hdelta, delta, file_ROC_lambda)
	hdelta = abs(hdelta);
	
	unique_value = unique(hdelta);
	unique_value = transpose(unique_value);

	tmp_result = [];
	for id = unique_value
		tp = nnz(hdelta > id & delta ~= 0);
		fp = nnz(hdelta > id & delta == 0);
		tn = nnz(hdelta <= id & delta == 0);
		fn = nnz(hdelta <= id & delta ~= 0);
		tmp_result = [tmp_result; [id, tp, fp, tn, fn]];
	end
	dlmwrite(file_ROC_lambda, tmp_result);
end
