function [P] = comp_P(N)
% generate P
% return a N*N matrix

	is = ([1:N]'-0.5)/N;

	js = zeros(1,N);

	Xs = is.^3/3-1/12;

	P = Xs+js;


end