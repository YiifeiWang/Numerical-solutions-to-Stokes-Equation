function [B] = spp_B(N)
% compute B in the saddle point problem
% return a 2N(N-1)*N^2 matrix

	B1 = comp_B1(N);
	B2 = comp_B2q(N);
	B = N*[B1;B2];

end