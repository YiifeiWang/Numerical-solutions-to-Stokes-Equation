addpath components

ns = [9,10,11];
vs = [20,40,80,160];
Ls = [2,4,6,8];

iter = zeros(3,5,4);
time = zeros(3,5,4);
res = zeros(3,5,4);
err = zeros(3,5,4);

for i = 1:length(ns)
	for j = 1+floor((i-1)/2):length(vs)
		for k = 1:length(Ls)
			if Ls(k)<ns(i)
				opts=[];
				opts.l = Ls(k);
				opts.max_num = 100;
				opts.tol = 1e-6;
				opts.smoother = 2;
				opts.v1 = vs(j);
				opts.v2 = vs(j);
				n = ns(i);
				out = VCycle_mod(n,opts);
				iter(i,j,k) = out.iter;
				time(i,j,k) = out.time;
				res(i,j,k) = out.res;
				u = out.U; v = out.V;
				err(i,j,k) = norm([u;v]-spp_U(2^n))/2^n;

				save('result/DSG_mod_result_tt_3.mat','iter','time','err','res');
				fprintf('N=%d, v1=v2=%d, l=%d, iter=%d, time=%.2f, err=%.2e, res =%.2e\n ',2^ns(i),vs(j),2^Ls(k),iter(i,j,k),time(i,j,k),err(i,j,k),res(i,j,k))
			end
		end
	end
end