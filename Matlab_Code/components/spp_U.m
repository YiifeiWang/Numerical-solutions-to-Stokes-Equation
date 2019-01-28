function [U] = spp_U(N)
% compute U in the saddle point problem
% return a 2N(N-1)*1 vector

	u = comp_U(N);
	v = comp_V(N);
	vt = v';
	U = [u(:);vt(:)];

end