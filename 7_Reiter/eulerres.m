function res = eulerres(kp_grid, k, j, kp, R, wage)
% A Simpel Euler equation residual function.
    global Params;
    % Unknown : kp
    % Current Node : k 
    LHS = (R*k + wage*Params.y_grid(j) - kp)^(-Params.mu);
    kpp = [interp1(Params.k_grid, kp_grid(:,1), kp), ...
           interp1(Params.k_grid, kp_grid(:,2), kp)];
    

    RHS = R * Params.beta * ((Params.Pz(j,1)*((R*kp + wage*Params.y_grid(1) - kpp(1))^(-Params.mu)))+... 
                             (Params.Pz(j,2)*((R*kp + wage*Params.y_grid(2) - kpp(2))^(-Params.mu))) );
    res = LHS - RHS;
end
 