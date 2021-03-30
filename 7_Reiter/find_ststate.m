function [X_ss, Rt, waget] = find_ststate
% This function find the steady-state level of collocation points, and
% knot points on distribution.

    global Params;
    kp_grid = repmat(Params.k_grid, 1, Params.n_y);

    % Make a guess of Steady-State price level
    R = 1.03; wage = 1; tol = 1e-3; diff = 100; 
    K = 3;
    while diff > tol 
        for iter = 1:Params.n_iter
            kp_grid_old = kp_grid;
            % Update policy function.
            kp_grid = solve_pol(kp_grid_old, R, wage);

            error = max(abs(kp_grid(:) - kp_grid_old(:)));
            if (error < Params.kp_error_tol)
                break;
            end
        end
        kp_star = kp_grid;
        % Compute transition matrix from the policy function.
        Q_ss	= compute_trans(kp_star);
        % Find stationary distribution.
        dist_ss = find_dist_ss(Q_ss);
        Kt = Params.knotDistrK'*dist_ss(1:999) ...
           + Params.knotDistrK'*dist_ss(1000:1998);
        L = 0.9300;
        KL_ratio = Kt/L;
        Rt = intr(KL_ratio, exp(0));
        waget = cwage(KL_ratio, exp(0));
        diff = abs(Rt - R);
        % Update price
        R = 0.6*Rt + 0.4*R;
    end

% Return: Collocation points on saving function (15 * 2)
%         Collocation points on distribution (999 * 2)
%         log(zss) = 0 (zss = 1) (1)
%         Collocation points on saving function with expectation (15 * 2)
    X_ss = [kp_star(:,1); kp_star(:,2); dist_ss; 0; kp_star(:,1); kp_star(:,2)];
end
