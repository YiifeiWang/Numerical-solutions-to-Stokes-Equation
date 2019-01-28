function [A_base] = comp_A_base(N)
% generate the base for kronecker product for A_N
% return a N*N matrix

	S = comp_S(N-1);
	Z = sparse(zeros(1,N));
	S1 = [S;Z];
	S2 = [Z;S];
	A_base = -S1+S2;
	
end