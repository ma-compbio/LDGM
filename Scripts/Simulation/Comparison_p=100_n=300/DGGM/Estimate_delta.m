function result = Estimate_delta1(hsigma1, hsigma2, lambda_seq, delta, file_ROC)
	result = [];
	delta = delta(triu(true(size(delta)),1));
	for lambda = lambda_seq
		file_ROC_lambda = strcat(file_ROC, '_lambda=', num2str(lambda),'.txt');
		hdelta = differential_graph(hsigma1, hsigma2, lambda);
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
