function [W1] = comp_W1(N)
% generate W_N^{(1)}
% return a N*2N matrix

	I = sparse(eye(N));
	W1 = kron(I,[1,1]);
end