function [B1] = comp_B1(N)
% generate B_N^{(1)}
% return a N(N-1)*N^2 matrix

	IN = sparse(eye(N));
	S = comp_S(N-1);
	B1 = kron(IN,S);
	
end