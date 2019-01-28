function [B2] = comp_B2(N)
% generate B_N^{(2)}
% return a N(N-1)*N^2 matrix

	IN = eye(N);
	S = comp_S(N-1);
	B2 = kron(S,sparse(IN(1,:)));
	for i =2:N
		B2 = [B2;kron(S,sparse(IN(i,:)))];
	end
	
end