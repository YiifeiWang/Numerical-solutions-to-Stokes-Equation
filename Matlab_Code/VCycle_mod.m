function [out]=VCycle(n, opts)
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%
	% This prgram implements modified V-Cycle multigrid method
	% for the saddle point problem.
	%
	% Input:
	%		n --- the lattice scale N = 2^n
	%		opts --- option structure with fields
	%				  l --- scale of the most coarse grid 2^l
	%		   smoother --- 0: DGS parallel
	%						2: DGS seq via Matlab (notes)
	%				 v1 --- number of iteration of smoothing in Stage 1
	%				 v2 --- number of iteration of smoothing in Stage 2
	%				tol --- tolerance of the norm of the residual
	%				tau --- parameter for Inexact Uzawa
	%			max_num --- maximum number of V-cycle
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
	if ~isfield(opts, 'max_num'); opts.max_num = 1500; end
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
	tqdm = opts.tqdm;

	assert(n-l>=1, 'size error')

	% construct and store the matrix which will be needed in V-Cycle
	N = 2^n;
	F = spp_F(N); F0 = F(1:N*(N-1)); G0 = F(1+N*(N-1):2*N*(N-1));
	U0 = zeros(N*(N-1),1); V0 = zeros(N*(N-1),1); P0 = zeros(N^2,1);
	U_his = {}; V_his = {}; P_his = {};
	F_his = {}; G_his = {}; R_his = {};
	A_his = {}; B1_his = {}; B2_his = {}; 
	if smoother ==0
		BTB_his = {};
	end
	RU_his = {}; LU_his = {}; LP_his = {}; RP_his = {};

	for k=0:n-l
		Nk = 2^(n-k);
		if k>0
			[LUk, LPk] = vc_lift(Nk);
			LU_his{k+1} = LUk; LP_his{k+1} = LPk;
		end
		A_his{k+1} = spp_subA(Nk); B1_his{k+1} = spp_B1(Nk); B2_his{k+1} = spp_B2(Nk);
		if smoother ==0
			BTB_his{k+1} = spp_BTB(Nk);
		end
		if k<n-l
			RU_his{k+1} = vc_restrict(Nk);
			RP_his{k+1} = vc_restrict_P(Nk);
		end
	end

	% start of the V-Cycle
	k = 0;
	Nk = 2^n;

	switch smoother
		case 0
			[Uk,Vk,Pk] = DGS(U0,V0,P0,F0,G0,A_his{k+1},B1_his{k+1},B2_his{k+1},BTB_his{k+1},v1);
		case 2
			[Uk,Vk,Pk] = DGS_seqM(U0,V0,P0,F0,G0,A_his{k+1},B1_his{k+1},B2_his{k+1},v1,Nk);
	end

	U_his{1}=Uk;
	V_his{1}=Vk;
	P_his{1}=Pk;
	F_his{1}=F0;
	G_his{1}=G0;
	R_his{1}=zeros(N^2,1);
	resF = F0-A_his{1}*Uk-B1_his{1}*Pk;
	resG = G0-A_his{1}*Vk-B2_his{1}*Pk;
	resR = -B1_his{1}'*Uk-B2_his{1}'*Vk;

	iter_num = 0;
	while norm([resF;resG;resR])>N^2*tol && iter_num<max_num
		Fk = RU_his{1}*resF;
		Gk = RU_his{1}*resG;
		Rk = RP_his{1}*resR;
		k = 1;
		% Stage 1
		while k<n-l
			Nk = 2^(n-k);
			Uk = zeros(Nk*(Nk-1),1); Vk = zeros(Nk*(Nk-1),1); Pk = zeros(Nk^2,1);
			switch smoother
				case 0
					[Uk, Vk, Pk] = DGS_mod(Uk,Vk,Pk,Fk,Gk,A_his{k+1},B1_his{k+1},B2_his{k+1},BTB_his{k+1},v1,Rk);
				case 2
					[Uk, Vk, Pk] = DGS_seqM_mod(Uk,Vk,Pk,Fk,Gk,A_his{k+1},B1_his{k+1},B2_his{k+1},v1,Nk,Rk);
			end
			U_his{k+1}=Uk;
			V_his{k+1}=Vk;
			P_his{k+1}=Pk;
			F_his{k+1}=Fk;
			G_his{k+1}=Gk;
			R_his{k+1}=Rk;
			subresF = Fk-A_his{k+1}*Uk-B1_his{k+1}*Pk;
			subresG = Gk-A_his{k+1}*Vk-B2_his{k+1}*Pk;
			subresR = Rk-B1_his{k+1}'*Uk-B2_his{k+1}'*Vk;
			Fk = RU_his{k+1}*subresF;
			Gk = RU_his{k+1}*subresG;
			Rk = RP_his{k+1}*subresR;
			k = k+1;
		end
		Nk = 2^l;
		Uk = zeros(Nk*(Nk-1),1); Vk = zeros(Nk*(Nk-1),1); Pk = zeros(Nk^2,1);
		switch smoother
			case 0
				[Uk, Vk, Pk] = DGS_mod(Uk,Vk,Pk,Fk,Gk,A_his{k+1},B1_his{k+1},B2_his{k+1},BTB_his{k+1},v1,Rk);
			case 2
				[Uk, Vk, Pk] = DGS_seqM_mod(Uk,Vk,Pk,Fk,Gk,A_his{k+1},B1_his{k+1},B2_his{k+1},v1,Nk,Rk);
		end
		% Stage 2
		k = n-l;
		while k>0
			Uk = U_his{k}+LU_his{k+1}*Uk;
			Vk = V_his{k}+LU_his{k+1}*Vk;
			Pk = P_his{k}+LP_his{k+1}*Pk;
			Nk = 2^(n-k+1);
			Fk = F_his{k};
			Gk = G_his{k};
			Rk = R_his{k};
			switch smoother
				case 0
					[Uk, Vk, Pk] = DGS_mod(Uk,Vk,Pk,Fk,Gk,A_his{k},B1_his{k},B2_his{k},BTB_his{k},v2,Rk);
				case 2
					[Uk, Vk, Pk] = DGS_seqM_mod(Uk,Vk,Pk,Fk,Gk,A_his{k},B1_his{k},B2_his{k},v2,Nk,Rk);
			end
			k = k-1;
		end
		resF = F0-A_his{1}*Uk-B1_his{1}*Pk;
		resG = G0-A_his{1}*Vk-B2_his{1}*Pk;
		resR = -B1_his{1}'*Uk-B2_his{1}'*Vk;
		U_his{1} = Uk;
		V_his{1} = Vk;
		P_his{1} = Pk;
		iter_num = iter_num+1;
		% print the result
		if tqdm && mod(iter_num,10)==0
			timeiter = toc;
			fprintf('Iter: %d finished, res: %.2e, time ellapsed: %.2f\n',iter_num, norm([resF;resG;resR])/N^2, timeiter);
		end

	end

	resF = F0-A_his{1}*Uk-B1_his{1}*Pk;
	resG = G0-A_his{1}*Vk-B2_his{1}*Pk;
	resR = -B1_his{1}'*Uk-B2_his{1}'*Vk;

	time = toc;
	out = [];
	out.time = time;
	out.iter = iter_num;
	out.U = Uk;
	out.V = Vk;
	out.P = Pk;
	out.res = norm([resF;resG;resR])/N^2;
end