function [U,V,P]=DGS_seq(U0,V0,P0,F,G,A,B1,B2,v1,N)
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
	
	DLA = tril(A);
	UA = triu(A,1);

	for i = 1:v1
		U = U+DLA\(F-A*U-B1*P);
		V = V+DLA\(G-A*V-B2*P);
		[U,V,P] = DGS_InnerM(N,U,V,P);
	end

end