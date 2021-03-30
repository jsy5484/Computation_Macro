function wage = cwage(KL_ratio,z)
% This function finds wage using Aggregate Info
  global Params;
  wage = z*(1-Params.alpha)*(KL_ratio)^Params.alpha;
end