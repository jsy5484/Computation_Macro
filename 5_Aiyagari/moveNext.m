function C_next = moveNext(a_grid, edg_grid, C_tilde, prod_grid, N_ass, N_prod, R, w)
    C_next = zeros(N_ass, N_prod);

    for j = 1:2
        for i = 1:N_ass
            % using Budget Constraint when binds
            if a_grid(i) < edg_grid(1, j)
                C_next(i, j) = R*a_grid(i) + w*prod_grid(j);
            
            % else using interpolation    
            elseif a_grid(i) >= edg_grid(N_ass, j)
                
                C_next(i, j)= ((edg_grid(N_ass,j) - a_grid(i))/(edg_grid(N_ass,j)-edg_grid(N_ass-1,j)))*C_tilde(N_ass-1, j)+ ...
                              ((a_grid(i)-edg_grid(N_ass-1,j))/(edg_grid(N_ass,j)-edg_grid(N_ass-1,j)))*C_tilde(N_ass, j);  
            else
                [low_adj, upp_adj, low, upp] = search(a_grid(i), edg_grid(:, j), N_ass);
                
                C_next(i, j)= ((upp_adj-a_grid(i))/(upp_adj-low_adj))*C_tilde(low, j)+ ...
                              ((a_grid(i)-low_adj)/(upp_adj-low_adj))*C_tilde(upp, j);  
            end  
        end
    end

end
