function [u] = true_u(x,y)
% compute true value of u at (x,y)

	u = (1-cos(2*pi*x))*sin(2*pi*y);

end