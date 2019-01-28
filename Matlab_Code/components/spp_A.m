function [A] = spp_A(N)
% compute A in the saddle point problem
% return a 2N(N-1)*2N(N-1) matrix

	AN = comp_A(N);
	I2 = sparse(eye(2));
	A = N^2*kron(I2,AN);

end