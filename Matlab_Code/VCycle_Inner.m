function [U,V] = VCycle_Inner(U0, V0, F0, G0, A_his, RU_his, LU_his, n, l, v1, v2, max_numVC, tau);
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%
	% This prgram implements V-Cycle multigridmethod
	% as Inner loop for Inexact Uzawa Iteration.
	%
	%
	%
	% Author: Yifei Wang 01/2019
	%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	

	assert(n-l>=1, 'size error')

	% construct and store the matrix which will be needed in V-Cycle
	N = 2^n;
	U_his = {}; V_his = {};
	F_his = {}; G_his = {};

	% start of the V-Cycle
	k = 0;
	Nk = 2^n;

	[Uk, Vk] = GS(U0, V0, F0, G0, A_his{k+1}, v1);

	U_his{1}=Uk;
	V_his{1}=Vk;
	F_his{1}=F0;
	G_his{1}=G0;
	resF = F0-A_his{1}*Uk;
	resG = G0-A_his{1}*Vk;

	iter_num = 0;
	while norm([resF;resG])>N^2*tau && iter_num<max_numVC
		Fk = RU_his{1}*resF;
		Gk = RU_his{1}*resG;
		k = 1;
		% Stage 1
		while k<n-l
			Nk = 2^(n-k);
			Uk = zeros(Nk*(Nk-1),1); Vk = zeros(Nk*(Nk-1),1);
			[Uk, Vk] = GS(Uk, Vk, Fk, Gk, A_his{k+1}, v1);
			U_his{k+1}=Uk;
			V_his{k+1}=Vk;
			F_his{k+1}=Fk;
			G_his{k+1}=Gk;
			subresF = Fk-A_his{k+1}*Uk;
			subresG = Gk-A_his{k+1}*Vk;
			Fk = RU_his{k+1}*subresF;
			Gk = RU_his{k+1}*subresG;
			k = k+1;
		end
		Nk = 2^l;
		Uk = zeros(Nk*(Nk-1),1); Vk = zeros(Nk*(Nk-1),1);
		[Uk, Vk] = GS(Uk, Vk, Fk, Gk, A_his{k+1}, v1);

		% Stage 2
		k = n-l;
		while k>0
			Uk = U_his{k}+LU_his{k+1}*Uk;
			Vk = V_his{k}+LU_his{k+1}*Vk;
			Nk = 2^(n-k+1);
			Fk = F_his{k};
			Gk = G_his{k};
			[Uk, Vk] = GS(Uk, Vk, Fk, Gk, A_his{k}, v2);
			k = k-1;
		end
		resF = F0-A_his{1}*Uk;
		resG = G0-A_his{1}*Vk;
		U_his{1} = Uk;
		V_his{1} = Vk;
		iter_num = iter_num+1;
	end

	U = Uk; V = Vk;
end