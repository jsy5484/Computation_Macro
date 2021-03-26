function h_star = h_find(k, kp, z, psi, theta, grate, delta)
    int = @(hr) (psi/(1-theta))*((k^theta)*(exp(z)*hr)^(1-theta) - grate*kp ...
        + (1-delta)*k) + ((hr-1)*(k^theta)*(exp(z)^(1-theta))*(hr^(-theta)));
    h_star = fsolve(int, 0.2);
    if h_star >= 1
        h_star = 0.9999;
    end
end