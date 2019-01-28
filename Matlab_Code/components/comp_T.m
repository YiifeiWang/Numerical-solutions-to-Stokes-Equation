function [T] = comp_T(N)
% generate T_N 
% return a N*N matrix

	T = gallery('tridiag', N, -1, 2, -1);
	
end