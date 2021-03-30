function MU = expected(a_pol_col, a_pol_L, a_pol_H, param)
    mu = param(1); 
    R = param(2);
    w = param(3);
    chi = param(4);
    l_low = param(5);
    l_high = param(6);

    syms a
    MU =  (((R*a_pol_col + w*l_low - double(subs(a_pol_L, a, a_pol_col)))^-mu) ...
        + ((R*a_pol_col + w*l_high - double(subs(a_pol_H, a, a_pol_col)))^-mu) ...
        + (chi*(min(a_pol_col, 0)^2)));
            
end