function [A] = comp_A(N)
% generate A_N
% return a N(N-1)*N(N-1) matrix

	A_base = comp_A_base(N);
	I = sparse(eye(N-1));
	IN = sparse(eye(N));
	T = comp_T(N-1);
	A = kron(IN,T)+kron(A_base,I);
	
end