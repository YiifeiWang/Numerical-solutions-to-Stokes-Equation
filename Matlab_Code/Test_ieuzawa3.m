addpath components

ns = [10,11];
vs = [2,4,8,16];
taus = [1e-4, 1e-5, 1e-6, 1e-7];
Ls = [2,4,6,8];

iter = zeros(2,4,4,4);
time = zeros(2,4,4,4);
err = zeros(2,4,4,4);
res = zeros(2,4,4,4);

for i = 1:length(ns)
	for j = 1:length(vs)
		for k = 1:length(taus)
			for l = 1:length(Ls) 
				opts=[];
				opts.Alpha = 1;
				opts.l = Ls(l);
				opts.tau = taus(k);
				opts.tol = 1e-7;
				opts.max_num = 10;
				opts.smoother = 4;
				opts.v1 = vs(j);
				opts.v2 = vs(j);
				n = ns(i);
				out = VCycle(n,opts);
				iter(i,j,k,l) = out.iter;
				time(i,j,k,l) = out.time;
				res(i,j,k,l) = out.res;
				u = out.U; v = out.V;
				err(i,j,k,l) = norm([u;v]-spp_U(2^n))/2^n;

				save('result/ieuzawa_result_tt_3.mat','iter','time','res','err');
				fprintf('N=%d, v1=v2=%d, tau=%f, l=%d, iter=%d, time=%.2f, err=%.2e, res =%.2e\n ',2^ns(i),vs(j),taus(k),2^Ls(l),iter(i,j,k,l),time(i,j,k,l),err(i,j,k,l),res(i,j,k,l))
			end
		end
	end
end