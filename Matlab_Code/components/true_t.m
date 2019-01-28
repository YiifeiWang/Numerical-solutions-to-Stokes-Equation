function [t] = true_t(x)
% compute true value of t at (x,1)

	t = 2*pi*(1-cos(2*pi*x));

end