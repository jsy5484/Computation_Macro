function res = res_equation(X, Xlag)
% This function takes Xt and Xt_1 and return residual 

    global Params;
    res = zeros(Params.n_states,1);
    
    % Inputs: Xt, Xt-1
    k =      X(1:60);
    lambda = X(61:2058);
    logZ =   X(2059);
    Ek =     X(2060:2119);
    
    k_lag =      Xlag(1:60);
    lambda_lag = Xlag(61:2058);
    logZ_lag =   Xlag(2059);
    Ek_lag =     Xlag(2060:2119);

    % Compute density
    Q = compute_trans(reshape(k, Params.n_k, Params.n_y));
    lambda_new = Q'*lambda_lag;
    
    % Solve price first.
    cdf_temp = lambda_lag(1:999) + lambda_lag(1000:1998);
    K = Params.knotDistrK'*cdf_temp;
    L = 0.5*(1-0.04)*1 + 0.5*(1-0.1)*1;
    KL_ratio = K/L;
    R_lag    = intr (KL_ratio, exp(logZ_lag));
    wage_lag = cwage(KL_ratio, exp(logZ_lag));
    
    clear cdf_temp K L KL_ratio
    cdf_temp = lambda_new(1:999) + lambda_new(1000:1998);
    K = Params.knotDistrK'*cdf_temp;
    L = 0.5*(1-0.04)*1 + 0.5*(1-0.1)*1;
    KL_ratio = K/L;
    R    = intr (KL_ratio, exp(Params.rho*logZ_lag));
    wage = cwage(KL_ratio, exp(Params.rho*logZ_lag));
    
    % Compute capital policy functions
    kp_grid_new = solve_pol_gen(reshape(Ek_lag, Params.n_k, Params.n_y), R_lag, wage_lag, R, wage);
    kp_grid_new = reshape(kp_grid_new, Params.n_kp, 1);

    
    % Compute residual
    res(1:60)      = k - kp_grid_new;
    res(2060:2119) = k - Ek_lag;
    res(61:2058)   = lambda - Q'*lambda_lag;
    res(2059)      = logZ - Params.rho*logZ_lag;
end

