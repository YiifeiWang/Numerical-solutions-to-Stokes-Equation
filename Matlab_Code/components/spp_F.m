function [F] = spp_F(N)
% compute F in the saddle point problem
% return a 2N(N-1)*1 vector

	f = comp_F(N);
	g = comp_G(N); gg = g';
	b = bnd_b(N); t = bnd_t(N);
	l = bnd_l(N); r = bnd_r(N);
	I = eye(N);
	I1 = I(:,1); IN = I(:,N);

	fbt = f(:)+N*(kron(I1,b)+kron(IN,t));
	glr = gg(:)+N*(kron(I1,l)+kron(IN,r));
	F = [fbt;glr];

end