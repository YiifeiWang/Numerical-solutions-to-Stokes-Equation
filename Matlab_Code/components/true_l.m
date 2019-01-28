function [l] = true_l(y)
% compute true value of l at (0,y)

	l = 2*pi*(1-cos(2*pi*y));

end