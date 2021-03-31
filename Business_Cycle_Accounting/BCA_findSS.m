function SS = BCA_findSS(Sbar)
global Params;
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

SS = [ks cs ls xs ys zs tauls tauxs gs];

end
