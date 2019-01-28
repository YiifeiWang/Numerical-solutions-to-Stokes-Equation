function [g] = true_g(x,y)
% compute true value of g at (x,y)

	g = 4*pi^2*(2*cos(2*pi*y)-1)*sin(2*pi*x);

end