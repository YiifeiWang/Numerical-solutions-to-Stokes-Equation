function [A] = spp_subA(N)
% compute A in the saddle point problem
% return a N(N-1)*N(N-1) matrix

	A = N^2*comp_A(N);

end