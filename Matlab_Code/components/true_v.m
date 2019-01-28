function [v] = true_v(x,y)
% compute true value of v at (x,y)

	v = -(1-cos(2*pi*y))*sin(2*pi*x);

end