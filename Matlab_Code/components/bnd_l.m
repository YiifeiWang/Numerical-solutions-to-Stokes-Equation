function [l] = bnd_l(N)
% generate l
% return a N*1 vector

	is = [1:N-1]'/N;
	l = 2*pi*(1-cos(2*pi*is));

end