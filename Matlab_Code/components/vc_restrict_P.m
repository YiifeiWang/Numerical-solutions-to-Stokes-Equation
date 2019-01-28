function [RP] = vc_restric_P(N)
% compute the restriction operator for U and F
% return a N^2/4*N^2 matrix

W1 = comp_W1(N/2);

RP = 1/4*kron(W1,W1);

end