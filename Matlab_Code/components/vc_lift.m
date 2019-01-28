function [LU, LP] = vc_lift(N)
% compute the lifting operator for U and P
% return a 4N(2N-1)*2N(N-1) matrix and a 2N(2N-1)*N(N-1)

W3 = comp_W1(N)';
W4 = comp_W2(N)';

LU = 1/2*kron(W3,W4);
%LU = kron(sparse(eye(2)),LU);

LP = kron(W3,W3);

end