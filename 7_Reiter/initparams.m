% Initialize parameters.
Params.beta  = 0.9;
Params.alpha = 0.33;
Params.delta = 0.1;
Params.rho   = 0.7;
Params.mu    = 2.0;
Params.sigma = 0.01;
Params.kmin  = 0.01; 
Params.kmax  = 30;
Params.ymax  = 1.01; 
Params.ymin  = 0.99;

Params.n_k    = 30;
Params.ndstst = 1000;
Params.n_y    = 2;            
Params.y_grid = [0.0, 1.0];
Params.Pz     = [0.5562, 0.4438; 0.0334, 0.9666];
Params.Qz     = kron(Params.Pz, ones(999,1));
Params.n_kp   = Params.n_y * Params.n_k;

Params.n_iter	     = 100;
Params.n_dist_iter	 = 1000;
Params.kp_error_tol	 = 1e-10;
Params.dist_error_tol= 1e-8;
Params.epsilon       = 1e-3;
Params.fsolve_options  = optimoptions('fsolve','Display','none','optimalityTolerance', 1e-15);
error_tol              = 1e-10;

