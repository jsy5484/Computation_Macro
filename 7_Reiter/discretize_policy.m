function Q = discretize_policy(grid, pol)
% Discretize policy function over the finer grid over distribution to construct Q

n_grid = length(grid);
    Q = zeros(length(pol),length(grid));

    for i = 1:length(pol)
        l_ind = find(pol(i) > grid, 1, 'last');
        if (isempty(l_ind))
            Q(i,1) = 1;
        elseif (l_ind == n_grid)
            Q(i,n_grid) = 1;
        else
            Q(i,l_ind)   = (grid(l_ind+1) - pol(i))/(grid(l_ind+1) - grid(l_ind));
            Q(i,l_ind+1) = 1 - Q(i,l_ind);
        end
    end
end
