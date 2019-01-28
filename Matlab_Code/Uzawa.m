function [U,V,P]=Uzawa(U0,V0,P0,F,G,A,B1,B2,v1,N,Alpha)
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%
	% This prgram implements Uzawa Iteration
	% for the saddle point problem.
	%
	%
	% Author: Yifei Wang 01/2019
	%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	U = U0; V = V0; P = P0;
	%FF = N^2*F; GG = N^2*G; AA = N^2*A; BB1 = N^2*B1; BB2 = N^2*B2;

	for i = 1:v1
		F1 = F-B1*P; G1 = G-B2*P;
		U = A\F1;
		V = A\G1;
		P = P+Alpha*(B1'*U+B2'*V);
	end

end