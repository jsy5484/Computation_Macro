function Q = compute_trans(kp_star)
% Compute an transition matrix Q with respect to all possible states.
    global Params;
    g_spl = [spline(Params.k_grid, kp_star(:,1), Params.knotDistrK),...
             spline(Params.k_grid, kp_star(:,2), Params.knotDistrK)];

    Qk = discretize_policy(Params.knotDistrK, reshape(g_spl, 999*2, 1));
    %Qk = discretize_policy(Params.knotDistrK, reshape(kp_star, 999*2, 1));

    Q = repelem(Params.Qz, 1, 999).*repmat(Qk, 1, 2);
end