function [x] = CG_solver(A,b,x0,tol)
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%
	% This prgram implements Conjudate Gradient Solver
	%
	%
	% Author: Yifei Wang 01/2019
	%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	if nargin<3; x0 = b; end;
	if nargin<4; tol = 1e-13; end;

	[N,~] = size(x0);
	r = b-A*x0;
	k = 0;
	x = x0;
	norm_rs = r'*r;
	alpha0 = norm_rs;
	tol2 = tol^2;
	while norm_rs>tol2;
		k = k+1;
		if k==1
			p=r;
		else
			betak = norm_rs/norm_rs_his;
			p = r + betak*p;
		end
		Ap = A*p;
		alphak = norm_rs/(p'*Ap);
		norm_rs_his = norm_rs;
		x = x + alphak*p;
		%r = b - A*x;
		r = r - alphak*Ap;
		norm_rs = r'*r;
	end
end