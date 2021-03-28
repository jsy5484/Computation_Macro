function [k_pol, c] = solveIND(Trans, unemp_gd, unemp_bd, ngrid_k, ngrid_KM, k_grid, KM_grid, k_pol, param, B, KM_min, KM_max, k_min, k_max)
    alpha = param(1);
    beta = param(2);
    delta = param(3);
    mu = param(4);
    
    pr_bd_unemp = zeros(ngrid_k, ngrid_KM, 2, 2);
    pr_bd_emp = zeros(ngrid_k, ngrid_KM, 2, 2);
    pr_gd_unemp = zeros(ngrid_k, ngrid_KM, 2, 2);
    pr_gd_emp = zeros(ngrid_k, ngrid_KM, 2, 2);

    pr_bd_unemp(:,:,1,1) = Trans(1,1)*ones(ngrid_k, ngrid_KM);
    pr_bd_unemp(:,:,1,2) = Trans(2,1)*ones(ngrid_k, ngrid_KM);
    pr_bd_unemp(:,:,2,1) = Trans(3,1)*ones(ngrid_k, ngrid_KM);
    pr_bd_unemp(:,:,2,2) = Trans(4,1)*ones(ngrid_k, ngrid_KM);
    
    pr_bd_emp(:,:,1,1) = Trans(1,2)*ones(ngrid_k, ngrid_KM);
    pr_bd_emp(:,:,1,2) = Trans(2,2)*ones(ngrid_k, ngrid_KM);
    pr_bd_emp(:,:,2,1) = Trans(3,2)*ones(ngrid_k, ngrid_KM);
    pr_bd_emp(:,:,2,2) = Trans(4,2)*ones(ngrid_k, ngrid_KM);
    
    pr_gd_unemp(:,:,1,1) = Trans(1,3)*ones(ngrid_k, ngrid_KM);
    pr_gd_unemp(:,:,1,2) = Trans(2,3)*ones(ngrid_k, ngrid_KM);
    pr_gd_unemp(:,:,2,1) = Trans(3,3)*ones(ngrid_k, ngrid_KM);
    pr_gd_unemp(:,:,2,2) = Trans(4,3)*ones(ngrid_k, ngrid_KM);
    
    pr_gd_emp(:,:,1,1) = Trans(1,4)*ones(ngrid_k, ngrid_KM);
    pr_gd_emp(:,:,1,2) = Trans(2,4)*ones(ngrid_k, ngrid_KM);
    pr_gd_emp(:,:,2,1) = Trans(3,4)*ones(ngrid_k, ngrid_KM);
    pr_gd_emp(:,:,2,2) = Trans(4,4)*ones(ngrid_k, ngrid_KM);
 
%% Calculate Price with aggregate labour and capital
    k_vec = zeros(ngrid_k, ngrid_KM, 2, 2);
    k_vec(:,:,1,1) = k_grid*ones(1, ngrid_KM);
    k_vec(:,:,1,2) = k_grid*ones(1, ngrid_KM);
    k_vec(:,:,2,1) = k_grid*ones(1, ngrid_KM);
    k_vec(:,:,2,2) = k_grid*ones(1, ngrid_KM);

    KM_vec = zeros(ngrid_k, ngrid_KM, 2, 2);
    KM_vec(:,:,1,1) = ones(ngrid_k, 1)*KM_grid';
    KM_vec(:,:,1,2) = ones(ngrid_k, 1)*KM_grid';
    KM_vec(:,:,2,1) = ones(ngrid_k, 1)*KM_grid';
    KM_vec(:,:,2,2) = ones(ngrid_k, 1)*KM_grid';
    
    agg_lab = zeros(ngrid_k, ngrid_KM, 2, 2);
    agg_lab(:,:,1,:) = (1-unemp_bd)*ones(ngrid_k, ngrid_KM, 2);
    agg_lab(:,:,2,:) = (1-unemp_gd)*ones(ngrid_k, ngrid_KM, 2);
    
    agg_shock_vec = zeros(ngrid_k, ngrid_KM, 2, 2);
    agg_shock_vec(:,:,1,:) = 0.99*ones(ngrid_k, ngrid_KM, 2);
    agg_shock_vec(:,:,2,:) = 1.01*ones(ngrid_k, ngrid_KM, 2);
    
    idio_shock_vec = zeros(ngrid_k, ngrid_KM, 2, 2);
    idio_shock_vec(:,:,:,1) = 0*ones(ngrid_k, ngrid_KM, 2);
    idio_shock_vec(:,:,:,2) = 1*ones(ngrid_k, ngrid_KM, 2);

    % Price from F.O.C. of Firm's problem
    one_mat = ones(ngrid_k, ngrid_KM, 2, 2);
    r_vec = alpha*(agg_shock_vec.*(KM_vec./agg_lab).^(alpha-1)) + (1-delta);
    w_vec = (1-alpha)*(agg_shock_vec.*(KM_vec./agg_lab).^alpha);
     
    wlth_vec = r_vec.*k_vec + w_vec.*idio_shock_vec;
    
    
%% Update the aggregate law of motion.
    KM_tmr = zeros(ngrid_k, ngrid_KM, 2, 2);
    KM_tmr(:,:,1,1) = exp(B(1) + B(2)*log(KM_vec(:,:,1,1)));
    KM_tmr(:,:,1,2) = exp(B(1) + B(2)*log(KM_vec(:,:,1,2)));
    KM_tmr(:,:,2,1) = exp(B(1) + B(2)*log(KM_vec(:,:,2,1)));
    KM_tmr(:,:,2,2) = exp(B(1) + B(2)*log(KM_vec(:,:,2,2)));
    KM_tmr = (KM_tmr >= KM_min).*(KM_tmr <= KM_max).*KM_tmr + ...
             (KM_tmr < KM_min)*KM_min + (KM_tmr > KM_max)*KM_max;
    
%% Construct path of future interest rate and wage
    r_b = 0.99*alpha*(KM_tmr./((1-unemp_bd)*one_mat)).^(alpha-1) + (1-delta);
    r_g = 1.01*alpha*(KM_tmr./((1-unemp_gd)*one_mat)).^(alpha-1) + (1-delta);
    
    w_b = 0.99*alpha*(KM_tmr./((1-unemp_bd)*one_mat)).^(alpha-1);
    w_g = 1.01*alpha*(KM_tmr./((1-unemp_gd)*one_mat)).^(alpha-1);

%% Now, solve our main problem
    diff = 100;
    tol = 1e-8;
    iter = 0;

    while diff>tol
        iter = iter + 1;
       %fprintf('Current Iteration is : %d \n', iter);

        % Recession & Unemployed
        kpp_bd_unemp = interpn(k_grid, KM_grid, k_pol(:,:,1,1), k_pol, KM_tmr, 'spline');
        consm_bd_unemp = r_b.*k_pol  - kpp_bd_unemp;
        consm_bd_unemp = (consm_bd_unemp>0).*consm_bd_unemp + (consm_bd_unemp<=0)*10^-10;
        MU_bd_unemp = consm_bd_unemp.^(-mu);
        
        % Recession & Employed
        kpp_bd_emp = interpn(k_grid, KM_grid, k_pol(:,:,1,2), k_pol, KM_tmr, 'spline');
        consm_bd_emp = r_b.*k_pol + (w_b.*one_mat) - kpp_bd_emp;
        consm_bd_emp = (consm_bd_emp>0).*consm_bd_emp + (consm_bd_emp<=0)*10^-10;
        MU_bd_emp = consm_bd_emp.^(-mu);
        
        % Boom & Unemployed
        kpp_gd_unemp = interpn(k_grid, KM_grid, k_pol(:,:,2,1), k_pol, KM_tmr, 'spline');
        consm_gd_unemp = r_g.*k_pol - kpp_gd_unemp;
        consm_gd_unemp = (consm_gd_unemp>0).*consm_gd_unemp + (consm_gd_unemp<=0)*10^-10;
        MU_gd_unemp = consm_gd_unemp.^(-mu);
        
        % Boom & Employed
        kpp_gd_emp = interpn(k_grid, KM_grid, k_pol(:,:,2,2), k_pol, KM_tmr, 'spline');
        consm_gd_emp = r_g.*k_pol  + (w_g.*one_mat) - kpp_gd_emp;
        consm_gd_emp = (consm_gd_emp>0).*consm_gd_emp + (consm_gd_emp<=0)*10^-10;
        MU_gd_emp = consm_gd_emp.^(-mu);
        
        % Expectation Operator
        RHS = (MU_bd_unemp.*r_b).*pr_bd_unemp + ...
              (MU_bd_emp.*r_b).*pr_bd_emp + ...
              (MU_gd_unemp.*r_g).*pr_gd_unemp +...
              (MU_gd_emp.*r_g).*pr_gd_emp;
        clear consm
        consm = (beta*RHS).^(-1/mu);
        k_pol_new = wlth_vec - consm;
        k_pol_new = (k_pol_new >= k_min).*(k_pol_new <= k_max).*k_pol_new + ...
                    (k_pol_new < k_min)*k_min + (k_pol_new > k_max)*k_max;
        
        diff = max(max(max(max(abs(k_pol_new - k_pol)))));
        k_pol = k_pol_new;
        
        %if diff<tol
        %    fprintf('Convergence obatined in: %d \n', iter);
        %end
    end
    
c = wlth_vec - k_pol;
