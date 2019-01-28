function [P] = spp_P(N)
% compute P in the saddle point problem
% return a N^2*1 vector

	p = comp_P(N);
	P = p(:);

end