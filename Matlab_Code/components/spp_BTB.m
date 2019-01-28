function [BTB] = ssp_BTB(N)
% generate B^TB for the saddle point problem
% return a N^2*N^2 matrix

	S2 = comp_S2(N);
	I = sparse(eye(N));
	BTB = N^2*(kron(S2,I)+kron(I,S2));

end