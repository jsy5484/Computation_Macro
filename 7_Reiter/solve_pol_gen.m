function kp_new = solve_pol_gen(kp_old, R_lag, wage_lag, R, wage)
global Params
kp_new = zeros(size(kp_old));

for i = 1:Params.n_k
	for j = 1:Params.n_y
        % Current nodes
		k = Params.k_grid(i);
		obj_fn = @(kp) eulerres_gen(kp_old, k, j, kp, R_lag, wage_lag, R, wage);

		% 	first check if lower bound binds
		diff_k = obj_fn(Params.k_grid(1));

		if (diff_k > 0)
			% lower bound binds
			kp_star = Params.kmin;
        else
            kp_star = fsolve(obj_fn, kp_old(i, j), Params.fsolve_options);
		end
		kp_new(i,j)	= kp_star;
	end
end
end