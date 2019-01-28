function [S2] = comp_S2(N)
% generate S_N^{(2)}
% return a N*(N+1) matrix

	S2 = gallery('tridiag', N, -1, 2, -1);
	S2(1,1) = 1;
	S2(N,N) = 1;

end