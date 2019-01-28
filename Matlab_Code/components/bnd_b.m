function [b] = bnd_b(N)
% generate b
% return a N*1 vector

	is = [1:N-1]'/N;
	b = -2*pi*(1-cos(2*pi*is));

end