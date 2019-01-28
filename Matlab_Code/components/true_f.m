function [f] = true_f(x,y)
% compute true value of f at (x,y)

	f = -4*pi^2*(2*cos(2*pi*x)-1)*sin(2*pi*y)+x^2;

end