function [extr_St, lkhat] = BCA_extract(C, Sbar, Yhat)
global Params;
z        = exp(Sbar(1)); 
taul     = Sbar(2);
taux     = Sbar(3);
g        = exp(Sbar(4));

kl       = ((1+taux)*(1-Params.beth*(1-Params.delta))/(Params.beth*Params.theta))^(1/(Params.theta-1))*z;
yk       = (kl/z)^(Params.theta-1);
xi1      = yk-(1+Params.gz)*(1+Params.gn)+1-Params.delta;
xi2      = (1-taul)*(1-Params.theta)*(kl)^Params.theta*z^(1-Params.theta)/Params.psi;
k        = (xi2+g)/(xi1+xi2/kl);
c        = xi1*k-g;
l        = k/kl;
y        = yk*k;
x        = y-c-g;

% log(ss)
lkss       = log(k);
lcss       = log(c);
lhss       = log(l);
lyss       = log(y);
lxss       = log(x);
lgss       = log(g);
lzss       = log(z);

yhat  = Yhat(:,1); ihat  = Yhat(:,2); hhat  = Yhat(:,3); ghat  = Yhat(:,4);
lyhat = log(yhat); lxhat = log(ihat); lhhat = log(hhat); lghat = log(ghat);
horiz = length(yhat);

% lnxt     = Y(:,6);
lkhat(1,1) = lkss;
Kthat(1,1) = exp(lkss);
 
for i=1:horiz
  lktp(i,1)    = lkss + ((1 - Params.delta)*(lkhat(i) - lkss) + x/k*(lxhat(i)-lxss))/(1+Params.gz)/(1+Params.gn);
  lkhat(i+1,1) = lktp(i);
  Ktp(i,1)     = ((1-Params.delta)*Kthat(i)+exp(lxhat(i)))/(1+Params.gz)/(1+Params.gn);
  Kthat(i+1,1) = Ktp(i);
end

lkhat    = lkhat(1:horiz);
Kthat    = Kthat(1:horiz);
lct      = lcss +(y*(lyhat-lyss) -x*(lxhat-lxss) -g*(lghat-lgss))/c;
lzt      = lzss +(lyhat -lyss -Params.theta*(lkhat-lkss))/(1 -Params.theta) -lhhat +lhss;
tault    = taul +(1-taul)*(lyhat-lyss-lct+lcss-1/(1-l)*(lhhat-lhss));
tauxt    = (lxhat-C(2,1)*lkhat-C(2,2)*lzt-C(2,3)*tault-C(2,5)*lghat-C(2,6))/C(2,4);

Ct       = exp(lyhat)-exp(lxhat)-exp(lghat);
Zt       = (exp(lyhat)./(Kthat.^Params.theta.*exp(lhhat).^(1-Params.theta))).^(1/(1-Params.theta));
Tault    = 1-Params.psi/(1-Params.theta)* (Ct./exp(lyhat)) .*(exp(lhhat)./(1-exp(lhhat)));
extr_St = [lzt, tault, tauxt, lghat, Tault, Zt];
end
