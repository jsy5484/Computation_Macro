function Val = Val_func(k) 
    global Params

    % do the interpolation
    V_interp = interp1(Params.k_grid, Params.V_old, k, 'linear');
    c = Params.z_now*Params.k_now^Params.theta + (1 - Params.delta)*Params.k_now - k; 
    
    if c <= 0
        Val = -9999999 - 999*abs(c);
    else
        Val = ((c^(1-Params.alpha))/(1 - Params.alpha)) + Params.beta*V_interp*Params.trans(:,Params.temp);
    end
        
    Val = - Val; 
    
end
