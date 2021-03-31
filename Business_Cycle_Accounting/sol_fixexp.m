function [A,B,C,X0] = sol_fixexp(Sbar,P,Q,s0,As)
% This function solves equilibrium of model under the control of expecation
% given Steady-state level, P, Q and P0. 
% Returns : A, B, and C matrices.

% Note : I followed Professor McGrattan's code very closedly


global Params;
P0      = (eye(4)-P)*Sbar;
param   = [Params.gn;Params.gz;Params.beta;Params.delta;Params.psi;Params.sigma;Params.theta;s0;As];

zs         = exp(Sbar(1));
tauls      = Sbar(2);
tauxs      = Sbar(3);
gs         = exp(Sbar(4));
kls        = ((1+tauxs)*(1-Params.beth*(1-Params.delta))/(Params.beth*Params.theta))^(1/(Params.theta-1))*zs;
A          = (zs/kls)^(1-Params.theta)-(1+Params.gz)*(1+Params.gn)+1-Params.delta;
B          = (1-tauls)*(1-Params.theta)*kls^Params.theta*zs^(1-Params.theta)/Params.psi;
ks         = (B+gs)/(A+B/kls);
cs         = A*ks-gs;
ls         = ks/kls;
ys         = ks^Params.theta*(zs*ls)^(1-Params.theta);
xs         = ys-cs-gs;
X0         = [log(ks);log(zs);tauls;tauxs;log(gs);1];
Y0         = [log(ys);log(xs);log(ls);log(gs)];

Z          = [log(ks);log(ks);log(ks);log(zs);log(zs);tauls;tauls;
              tauxs;tauxs;log(gs);log(gs)];
del        = max(abs(Z)*1e-5,1e-8);
for i=1:11;
  Zp       = Z;
  Zm       = Z;
  Zp(i)    = Z(i)+del(i);
  Zm(i)    = Z(i)-del(i);
  dR(i,1)  = (eulerres2(Zp,param)-eulerres2(Zm,param))/(2*del(i));
end;

a0         = dR(1);
a1         = dR(2);
a2         = dR(3);
b0         = dR(4:2:11)';
b1         = dR(5:2:11)';
tem        = roots([a0,a1,a2]);
gammak     = tem(find(abs(tem)<1));
gamma      = -((a0*gammak+a1)*eye(4)+a0*P')\(b0*P+b1)';
gamma0     = (1-gammak)*log(ks)-gamma'*[log(zs);tauls;tauxs;log(gs)];
Gamma      = [gammak;gamma;gamma0];

philh      =-(Params.psi*ys*(1-Params.theta)+(1-Params.theta)*(1-tauls)*ys*(1-ls)/ls*Params.theta+ ...
             (1-Params.theta)*(1-tauls)*ys);
philk      = (Params.psi*ys*Params.theta+Params.psi*(1-Params.delta)*ks- ...
             (1-Params.theta)*(1-tauls)*ys*(1-ls)/ls *Params.theta)/philh;
philz      = (Params.psi*ys*(1-Params.theta)-(1-Params.theta)^2*(1-tauls)*ys*(1-ls)/ls)/philh;
phill      = ((1-Params.theta)*(1-tauls)*ys*(1-ls)/ls *(1/(1-tauls)))/philh;
philg      = (-Params.psi*gs)/philh;
philkp     = (-Params.psi*(1+Params.gz)*(1+Params.gn)*ks)/philh;
phiyk      = Params.theta+(1-Params.theta)*philk;
phiyz      = (1-Params.theta)*(1+philz);
phiyl      = (1-Params.theta)*phill;
phiyg      = (1-Params.theta)*philg;
phiykp     = (1-Params.theta)*philkp;
phixk      = -ks/xs*(1-Params.delta);
phixkp     = ks/xs*(1+Params.gz)*(1+Params.gn);

A     = [ gammak,                         gamma', gamma0;
         [0;0;0;0],            P,                     P0;
               0,      0,      0,      0,      0,      1];
B     = [ 0,0,0,0;
            Q;
          0,0,0,0];
C     = [ [phiyk,phiyz*As(1),phiyl*As(2),    0,phiyg*As(4)]+phiykp*Gamma(1:5)';
          [phixk,       0,       0,    0,                0]+phixkp*Gamma(1:5)';
          [philk,philz*As(1),phill*As(2),    0,philg*As(4)]+philkp*Gamma(1:5)';
          0,0,0,0,As(4)];
phi0  = Y0-C*X0(1:5);
C     = [C,phi0];

