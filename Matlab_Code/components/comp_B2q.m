function [B2] = comp_B2(N)
% generate B_N^{(2)}
% return a N(N-1)*N^2 matrix
	
	aux1 = [1:(N-1)*N]';
	iindex = [aux1;aux1];
	aux2 = [0:(N-2)]'*N+[1:N];
	jindex = [aux2(:);aux2(:)+N];
	aux3 = zeros(N*(N-1),1);
	value = [aux3-1;aux3+1];
	B2 = sparse(iindex,jindex,value,N*(N-1),N^2);
	
end