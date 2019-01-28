function [W2] = comp_W2(N)
% generate W_N^{(2)}
% return a (N-1)*(2N-1) matrix

	W1 = comp_W1(N-1);
	Z = sparse(zeros(N-1,1));
	W2 = [W1,Z]+[Z,W1];
end