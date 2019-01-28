function [V] = comp_V(N)
% generate V
% return a (N-1)*N matrix

	is = ([1:N]-0.5)'/N;

	js = [1:N-1]/N;

	Ys = 1-cos(2*pi*js);
	Xs = sin(2*pi*is);

	V = -Xs.*Ys;


end