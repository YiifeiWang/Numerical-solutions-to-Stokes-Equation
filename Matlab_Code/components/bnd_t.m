function [t] = bnd_t(N)
% generate t
% return a N*1 vector

	is = [1:N-1]'/N;
	t = 2*pi*(1-cos(2*pi*is));

end