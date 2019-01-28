function [B] = spp_B(N)
% compute B1 in the saddle point problem
% return a N(N-1)*N^2 matrix

	B1 = comp_B1(N);
	B = N*B1;

end