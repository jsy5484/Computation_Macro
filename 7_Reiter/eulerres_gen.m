function res = eulerres_gen(kp_grid, k, j, kp, R_lag, wage_lag, R, wage)
    global Params;
    % Unknown : kp
    % Current Node : k 
    LHS = (R_lag*k + wage_lag*Params.y_grid(j) - kp)^(-Params.mu);
    kpp = [interp1(Params.k_grid, kp_grid(:,1), kp), ...
           interp1(Params.k_grid, kp_grid(:,2), kp)];
    

    RHS = R_lag * Params.beta * ((Params.Pz(j,1)*((R*kp + wage*Params.y_grid(1) - kpp(1))^(-Params.mu)))+... 
                             (Params.Pz(j,2)*((R*kp + wage*Params.y_grid(2) - kpp(2))^(-Params.mu))) );
    res = LHS - RHS;
end
 