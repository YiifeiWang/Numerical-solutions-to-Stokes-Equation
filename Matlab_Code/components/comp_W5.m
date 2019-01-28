function [W5] = comp_W5(N)
% generate W_N^{(5)}
% return a 2N*N matrix
	I = sparse(eye(N-1));
	Z = sparse(2*N-2,1);
	Z2 = sparse(1,N);

	aux1 = kron(I,[3;1]);
	aux2 = kron(I,[1;3]);

	W5 = [Z2;[aux1,Z]+[Z,aux2];Z2];
	W5(1,1) = 3;
	W5(2*N,N) = 3;
	
end