function [lambda_min] = Find_Lambda_min(hsigma1, hsigma2, lambda_max, n_target)
	lambda_min = lambda_max;
	tp = 0;
	iter = 0;
	while tp < n_target & iter < 10
		lambda_min = 0.5 * lambda_min;
		[hdelta] = differential_graph(hsigma1,hsigma2,lambda_min);
		tp = nnz(triu(hdelta, 1));
		iter = iter + 1;
	end
end

