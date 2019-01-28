function [r] = bnd_r(N)
% generate r
% return a N*1 vector

	is = [1:N-1]'/N;
	r = -2*pi*(1-cos(2*pi*is));

end