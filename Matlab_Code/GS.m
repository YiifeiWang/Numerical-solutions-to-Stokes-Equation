function [U,V]=GS(U0,V0,F,G,A,v1)
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%
	% This prgram implements Gauss-Seidel Iteration
	% as Inner loop of V-Cycle multi-grid.
	%
	%
	% Author: Yifei Wang 01/2019
	%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	U = U0; V = V0;
	
	DLA = tril(A);
	UA = triu(A,1);

	for i = 1:v1
		U = U+DLA\(F-A*U);
		V = V+DLA\(G-A*V);
	end

end