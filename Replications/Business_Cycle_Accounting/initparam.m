
Params.theta    = 0.35;
Params.delta    = 1-(1-.0464)^(1/4); 
Params.gz       = (1.016)^(1/4)-1; 
Params.gn       = (1.015)^(1/4)-1;
Params.beta     = 0.9722^(1/4); 
Params.psi      = 2.24;
Params.grate    = (1+Params.gn)*(1+Params.gz); 
Params.sigma    = 1.000001000000000;
Params.beth     = Params.beta*(1+Params.gz)^(-Params.sigma);
Params.adjb     = (1+Params.gn)*(1+Params.gz)-1+Params.delta;
Params.qadj     = 0.25/sum([Params.gn; Params.gz; Params.delta]);

Params.theta = 0.3333334; Params.delta = 0.05; Params.gamma_z = 0.016; Params.gamma_n = 0.015;
Params.beta = 0.9722; beta_hat = beta*(1+gamma_n); Params.psi = 2.5;
Params.grate    = (1+Params.gn)*(1+Params.gz); 
Params.sigma    = 1.000001000000000;
Params.beth     = Params.beta*(1+Params.gz)^(-Params.sigma);
Params.adjb     = (1+Params.gn)*(1+Params.gz)-1+Params.delta;
Params.qadj     = 0.25/sum([Params.gn; Params.gz; Params.delta]);