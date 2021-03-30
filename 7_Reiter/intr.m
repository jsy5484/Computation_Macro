function R = intr(KL_ratio,z)
% This function finds interest rate using Aggregate Info
    global Params;
    R = 1+z*Params.alpha*(KL_ratio)^(Params.alpha-1) - Params.delta;
end
