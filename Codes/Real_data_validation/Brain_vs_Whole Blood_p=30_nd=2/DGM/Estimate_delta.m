function result = Estimate_delta1(hsigma1, hsigma2, lambda_seq, delta, delta_union)
	result = [];
	for lambda = lambda_seq
		hdelta = differential_graph(hsigma1, hsigma2, lambda);
		[x, y, z] = Compare_delta_hdelta(hdelta, delta, delta_union);
		result = [result; [x, y, x, z]];
	end
end

function [x, y, z] = Compare_delta_hdelta(hdelta, delta, delta_union)
	A = hdelta(triu(true(size(hdelta)),1));
	x = nnz(A);
	
	idx = sub2ind(size(hdelta), delta(:,1), delta(:,2));
	B = hdelta(idx);
	y = nnz(B);

	idx = sub2ind(size(hdelta), delta_union(:,1), delta_union(:,2));
	C = hdelta(idx);
	z = nnz(C);
end

