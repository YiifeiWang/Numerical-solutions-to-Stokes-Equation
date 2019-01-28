function [G] = comp_G(N)
% generate G
% return a N*(N-1) matrix

	is = ([1:N]-0.5)'/N;

	js = [1:N-1]/N;

	Ys = 2*cos(2*pi*js)-1;
	Xs = sin(2*pi*is);

	G = 4*pi^2*Xs.*Ys;


end