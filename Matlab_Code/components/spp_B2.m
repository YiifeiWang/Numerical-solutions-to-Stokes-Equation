function [B] = spp_B2(N)
% compute B2 in the saddle point problem
% return a N(N-1)*N^2 matrix

	B2 = comp_B2q(N);
	B = N*B2;

end