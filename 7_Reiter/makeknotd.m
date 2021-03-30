% Matlab file for_ paper "Solving Heterogeneous Agent Models by Projection and Perturbation"
% Michael Reiter, Institute for Advanced Studies, September 2006
% Last update: June 2008
% Feel free to use, copy and modify at your own risk;
%   this program comes with NO WARRANTY WHATSOEVER
function [knots,logshift] = makeknotd(kmin,kmax,n,logshift);
  if(nargin>3)
    knots = logspaceshift(kmin,kmax,n,logshift);
  else
    [knots,logshift] = logspaceshift(kmin,kmax,n,1,n/4);
  end
  knots = knots';
