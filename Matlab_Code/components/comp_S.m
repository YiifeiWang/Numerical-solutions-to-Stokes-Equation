function [S] = comp_S(N)
% generate S_N
% return a N*(N+1) matrix

	S1 = gallery('tridiag', N, 0, -1, 1);
	I = eye(N);
	S2 = sparse(I(:,N));
	S=[S1,S2];
end