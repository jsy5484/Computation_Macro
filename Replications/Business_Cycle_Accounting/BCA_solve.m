function [A, B, C] = BCA_solve(SS, param, P, Q, P0, X0, Y0)
% This function solves equilibrium of model given Steady-state level, P, Q
% and P0. Returns : A, B, and C matrices.

% Note : I followed Professor McGrattan's code very closedly

global Params;
ks = SS(1); cs = SS(2); hs = SS(3); xs = SS(4); ys = SS(5); zs = SS(6);
tauls = SS(7); tauxs = SS(8); gs = SS(9);


Z          = [log(ks); log(ks); log(ks); log(zs); log(zs); tauls; tauls;
              tauxs; tauxs; log(gs); log(gs)];
del        = max(abs(Z)*1e-5, 1e-8);

% Calculate Numerical Derivatives around SS
for i=1:11
  Zp       = Z;
  Zm       = Z;
  Zp(i)    = Z(i) + del(i);
  Zm(i)    = Z(i) - del(i);
  dR(i,1)  = (eulerres(Zp, param) - eulerres(Zm, param))/(2*del(i));
end
 
a0         = dR(1);
a1         = dR(2);
a2         = dR(3);
b0         = dR(4:2:11)';
b1         = dR(5:2:11)';

temp        = roots([a0,a1,a2]);
gammak     = temp(find(abs(temp)<1));

if sum(size(gammak))~=2 || isreal(gammak) == 0
    L = 1e+20; dL = 1e+20*ones(size(Params.Theta));
    return
else
    gamma   = -((a0*gammak+a1)*eye(4) + a0*P')\(b0*P + b1)';
    gamma0  = (1 - gammak)*log(ks) -gamma'*[log(zs); tauls; tauxs; log(gs)];
    Gamma   = [gammak; gamma; gamma0];

    %
    %     State-space system:   X[t+1] = A X[t] + B eps[t+1]
    %                           Y[t]   = C X[t] + ome[t]
    %
    
    philh   = -(Params.psi*ys*(1 -Params.theta) + (1 -Params.theta)*(1 -tauls)*ys*(1-hs)/hs*Params.theta+ ...
                 (1 -Params.theta)*(1 -tauls)*ys);
    philk   = (Params.psi*ys*Params.theta +Params.psi*(1 -Params.delta)*ks- ...
                 (1 -Params.theta)*(1 -tauls)*ys*(1 -hs)/hs *Params.theta)/philh;
    philz   = (Params.psi*ys*(1 -Params.theta) -(1 -Params.theta)^2*(1 -tauls)*ys*(1 -hs)/hs)/philh;
    phill   = ((1 -Params.theta)*(1-tauls)*ys*(1-hs)/hs *(1/(1-tauls)))/philh;
    philg   = ( -Params.psi*gs)/philh;
    philkp  = ( -Params.psi*(1 +Params.gz)*(1 +Params.gn)*ks)/philh;
    phiyk   = Params.theta+(1 -Params.theta)*philk;
    phiyz   = (1 -Params.theta)*(1+philz);
    phiyl   = (1 -Params.theta)*phill;
    phiyg   = (1 -Params.theta)*philg;
    phiykp  = (1 -Params.theta)*philkp;
    phixk   = -ks/xs*(1 -Params.delta);
    phixkp  = ks/xs*(1 +Params.gz)*(1 +Params.gn);


    A     = [ gammak,                         gamma', gamma0;
             [0;0;0;0],            P,                     P0;
                   0,      0,      0,      0,      0,      1];
    B     = [ 0,0,0,0;
                Q;
              0,0,0,0];
    C     = [ [phiyk,phiyz,phiyl,    0,phiyg]+phiykp*Gamma(1:5)';
              [phixk,    0,    0,    0,    0]+phixkp*Gamma(1:5)';
              [philk,philz,phill,    0,philg]+philkp*Gamma(1:5)';
              0,0,0,0,1];
    phi0  = Y0-C*X0(1:5);
    C     = [C,phi0];
end