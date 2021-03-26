%% This is a small function finds a steady state of system

function [kss, hss] = SS(param)
    F = @(x) [ param(5)*(param(1)*(x(1)^(param(1)-1))*(x(2)^(1-param(1))) + 1-param(2))- (1+param(4))*(1+param(3));
              (-(param(6)/(1-x(2)))) + ((1-param(1))*(x(1)^param(1))*(x(2)^(-param(1)))) / ((x(1)^param(1))*(x(2)^(1-param(1))) + (1-param(2))*x(1) - (1+param(4))*(1+param(3))*x(1))];
    x0 = [0.5; 0.5];
    
    [kss, hss] = fsolve(F,x0);
end
