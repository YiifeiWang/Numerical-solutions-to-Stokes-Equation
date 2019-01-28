function [F] = comp_F(N)
% generate F
% return a (N-1)*N matrix

	is = [1:N-1]'/N;

	js = ([1:N]-0.5)/N;

	Xs = 2*cos(2*pi*is)-1;
	Ys = sin(2*pi*js);

	F = -4*pi^2*Xs.*Ys+is.^2;


end