function [lambda_max] = Find_Lambda_max(hsigma1, hsigma2)
	lambda_init = 1;
	[hdelta] = differential_graph(hsigma1,hsigma2,lambda_init);
	% ind = ~ any(hdelta(:));
	iter = 0;
	if nnz(triu(hdelta, 1)) == 0
		while nnz(triu(hdelta, 1)) == 0 & iter < 30
			lambda_init = lambda_init / 2;
			[hdelta] = differential_graph(hsigma1,hsigma2,lambda_init);
			iter = iter + 1;
		end
		lambda_max = lambda_init * 5;
	end

	if nnz(triu(hdelta, 1)) ~= 0
		while nnz(triu(hdelta, 1)) ~= 0 & iter < 30
			lambda_init = lambda_init * 2;
			[hdelta] = differential_graph(hsigma1,hsigma2,lambda_init);
			iter = iter + 1;
		end
		lambda_max = lambda_init;
	end
end
