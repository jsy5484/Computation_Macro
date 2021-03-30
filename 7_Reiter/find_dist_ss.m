function dist = find_dist_ss(transition_matrix)
    global Params;
    n_dist = 999 * Params.n_y;

    % Initial Guess 
    dist = 1./n_dist * ones(n_dist,1);

    for dist_iter = 1:Params.n_dist_iter
        dist_old = dist;
        dist = transition_matrix' * dist_old;

        error = max(abs(dist(:) - dist_old(:)));
        if (error < Params.dist_error_tol)
            break;
        end
    end
end