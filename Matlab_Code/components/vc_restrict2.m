function [RU] = vc_restric(N)
% compute the restriction operator for U and F
% return a N(N/2-1)*2N(N-1) matrix and a N/2(N/2-1)*N(N-1)

W1 = comp_W5(N/2)';
W2 = comp_W2(N/2);

RU = 1/32*kron(W1,W2);
%RU = kron(sparse(eye(2)),RU);

end