function [U,V,P]=DGS(U0,V0,P0,F,G,A,B1,B2,BTB,v1,R)
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%
	% This prgram implements Distributive Gauss-Seidel Iteration
	% for the saddle point problem.
	%
	%
	% Author: Yifei Wang 01/2019
	%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	U = U0; V = V0; P = P0;
	[~, N2] = size(B1);
	D = spdiags(1./diag(BTB),0,N2,N2);
	
	DLA = tril(A);
	UA = triu(A,1);

	for i = 1:v1
		U = U+DLA\(F-A*U-B1*P);
		V = V+DLA\(G-A*V-B2*P);
		Q = D*(R-B1'*U-B2'*V);
		U = U+B1*Q;
		V = V+B2*Q;
		P = P-BTB*Q;
	end

end