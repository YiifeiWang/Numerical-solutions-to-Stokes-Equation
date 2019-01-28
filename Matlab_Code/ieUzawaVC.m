function [out]=ieUzawa(n, opts)
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%
	% This prgram implements Inexact Uzawa Iteration 
	% for the saddle point problem based on V-Cycle 
	% multigrid method
	%
	% Input:
	%		n --- the lattice scale N = 2^n
	%		opts --- option structure with fields
	%				  l --- scale of the most coarse grid 2^l
	%				 v1 --- number of iteration of smoothing in Stage 1
	%				 v2 --- number of iteration of smoothing in Stage 2
	%				tol --- tolerance of the norm of the residual
	%				tau --- parameter for V-cycle
	%			max_num --- maximum number of outter iteration
	%		  max_numVC --- maximum number of V-cycle
	%			   tqdm --- whether to print the result
	%
	%
	% Author: Yifei Wang 01/2019
	%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	tic;
	if nargin < 2; opts = []; end
	if ~isfield(opts, 'smoother'); opts.smoother = 0; end
	if ~isfield(opts, 'l'); opts.l = 2; end
	if ~isfield(opts, 'v1'); opts.v1 = 100; end
	if ~isfield(opts, 'v2'); opts.v2 = 100; end
	if ~isfield(opts, 'Alpha'); opts.Alpha = 1; end
	if ~isfield(opts, 'tol'); opts.tol = 1e-8; end
	if ~isfield(opts, 'tau'); opts.tau = 1e-4; end
	if ~isfield(opts, 'max_num'); opts.max_num = 100; end
	if ~isfield(opts, 'max_numVC'); opts.max_numVC = 100; end
	if ~isfield(opts, 'tqdm'); opts.tqdm = 0; end

	% copy parameter
	smoother = opts.smoother;
	l = opts.l;
	v1 = opts.v1;
	v2 = opts.v2;
	Alpha = opts.Alpha;
	tol = opts.tol;
	tau = opts.tau;
	max_num = opts.max_num;
	max_numVC = opts.max_numVC;
	tqdm = opts.tqdm;

	assert(n-l>=1, 'size error')

	% construct and store the matrix which will be needed in V-Cycle
	N = 2^n;
	F = spp_F(N); F0 = F(1:N*(N-1)); G0 = F(1+N*(N-1):2*N*(N-1));
	U0 = zeros(N*(N-1),1); V0 = zeros(N*(N-1),1); P0 = zeros(N^2,1);
	A_his = {}; 
	RU_his = {}; LU_his = {};

	for k=0:n-l
		Nk = 2^(n-k);
		if k>0
			[LUk, LPk] = vc_lift(Nk);
			LU_his{k+1} = LUk;
		end
		A_his{k+1} = spp_subA(Nk);
		if k<n-l
			RU_his{k+1} = vc_restrict(Nk); 
		end
	end

	B1 = spp_B1(N); B2 = spp_B2(N);

	% start of the iteration
	outiter = 0;
	U = U0; V = V0; P = P0;
	F = F0 - B1*P; G = G0 - B2*P;
	[U,V] = VCycle_Inner(U, V, F, G, A_his, RU_his, LU_his, n, l, v1, v2, max_numVC, tau);
	P = P+Alpha*(B1'*U+B2'*V);

	resF = F0-A_his{1}*U-B1*P;
	resG = G0-A_his{1}*V-B2*P;

	outiter = outiter + 1;

	while norm([resF;resG])/N^2>tol && outiter<max_num
		F = F0 - B1*P; G = G0 - B2*P;
		[U,V] = VCycle_Inner(U, V, F, G, A_his, RU_his, LU_his, n, l, v1, v2, max_numVC, tau);
		P = P+Alpha*(B1'*U+B2'*V);

		resF = F0-A_his{1}*U-B1*P;
		resG = G0-A_his{1}*V-B2*P;

		outiter = outiter + 1;
		% print the result
		if tqdm && mod(outiter,10)==0
			timeiter = toc;
			fprintf('Iter: %d finished, res: %.2e, time ellapsed: %.2f\n',outiter, norm([resF;resG])/N^2, timeiter);
		end

	end

	U_his = {}; V_his = {};
	F_his = {}; G_his = {}; 

	time = toc;
	out = [];
	out.time = time;
	out.iter = outiter;
	out.U = U;
	out.V = V;
	out.P = P;
	out.res = norm([resF;resG])/N^2;
end