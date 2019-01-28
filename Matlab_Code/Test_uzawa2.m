addpath components

ns = [8,9];
vs = [2,4,8,16];
alphas = [0.5,0.75,1,1.5,2];
Ls = [2,4,6];

iter = zeros(2,4,4,3);
time = zeros(2,4,4,3);
err = zeros(2,4,4,3);
res = zeros(2,4,4,3);

for i = 1:length(ns)
	for j = 1:length(vs)
		for k = 1:length(alphas)
			for l = 1:length(Ls) 
				opts=[];
				opts.Alpha = alphas(k);
				opts.l = Ls(l);
				opts.max_num = 20;
				%opts.tol = 1e-7;
				opts.smoother = 3;
				opts.v1 = vs(j);
				opts.v2 = vs(j);
				n = ns(i);
				out = VCycle(n,opts);
				iter(i,j,k,l) = out.iter;
				time(i,j,k,l) = out.time;
				res(i,j,k,l) = out.res;
				u = out.U; v = out.V;
				err(i,j,k,l) = norm([u;v]-spp_U(2^n))/2^n;

				save('result/uzawa_result_tt_2.mat','iter','time','res','err');
				fprintf('N=%d, v1=v2=%d, alpha=%f, l=%d, iter=%d, time=%.2f, err=%.2e, res =%.2e\n ',2^ns(i),vs(j),alphas(k),2^Ls(l),iter(i,j,k,l),time(i,j,k,l),err(i,j,k,l),res(i,j,k,l))
			end
		end
	end
end