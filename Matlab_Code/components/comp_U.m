function [U] = comp_U(N)
% generate U
% return a (N-1)*N matrix

	is = [1:N-1]'/N;

	js = ([1:N]-0.5)/N;

	Xs = 1-cos(2*pi*is);
	Ys = sin(2*pi*js);

	U = Xs.*Ys;


end