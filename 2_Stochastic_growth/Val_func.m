function Val = Val_func(k) 
    global Params
    
    %kl_idx = max(sum(k>Params.kgrid), 1); % identify the gridpoint that falls just below . .
    %kh_idx = klo + 1;

    % do the interpolation
    %V_interp = V_old(kl_idx, :) + (k - kgrid(klo))*(v0(kh_idx, :) - v0(klo, :))/(kmat(kh_idx) - kmat(kl_idx));
    V_interp = interp1(Params.k_grid, Params.V_old, k, 'linear');
    
    c = Params.z_now*Params.k_now^Params.theta + (1 - Params.delta)*Params.k_now - k; 
    if c <= 0
        Val = -9999999 - 999*abs(c);
    else
        Val = ((c^(1-Params.alpha))/(1 - Params.alpha)) + Params.beta*V_interp*Params.trans(:,Params.temp);
    end
        
    Val = - Val; % make it negative since we're maximizing and code is to minimize.
    
end
