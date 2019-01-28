function [U,V,P]=Inexact_Uzawa(U0,V0,P0,F,G,A,B1,B2,v1,N,Alpha,tau)
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%
	% This prgram implements Inexact Uzawa Iteration
	% for the saddle point problem.
	%
	%
	% Author: Yifei Wang 01/2019
	%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	U = U0; V = V0; P = P0;

	for i = 1:v1
		F1 = F-B1*P; G1 = G-B2*P;
		U = CG_solver(A,F1,U,tau);
		V = CG_solver(A,G1,V,tau);
		P = P+Alpha*(B1'*U+B2'*V);
	end

end