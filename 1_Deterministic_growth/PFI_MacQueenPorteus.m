clear all clc
tic;

% Parameterization
alpha = 0.3333333333; delta = 0.1; beta = 0.975; A = 1.05;
maxDiff = 10.0; tol = 1.0e-10; iter = 0; nkgrid = 1001; m = 100;

% Find the Steady-state
kss   = ((1 + beta * (delta - 1))/(alpha * beta * A))^(1 / (alpha - 1));
yss   = A * kss^alpha + (1 - delta)*kss;
css   = yss - kss;

% Construct Capital Grid
kmin = 0.5*kss; kmax = 1.5*kss;
k_grid = linspace(kmin, kmax, nkgrid);
kp_grid = linspace(kmin, kmax, nkgrid);

kp_pol = k_grid;

%rtrn = zeros(ngrid_k, ngrid_k);
V_old = zeros(1, nkgrid);
V_nxt = zeros(1, nkgrid);
dr    = ones(1, nkgrid);

while maxDiff > tol
    iter = iter + 1;
    fprintf('Current Iteration is: %.d \n', iter);
    
    for j = 1:m
        for i = 1:nkgrid
            V_nxt(i) = log(A*(k_grid(i)^alpha) + (1 - delta)*k_grid(i) - kp_pol(i)) + beta*V_old(dr(i));
        end
        
        V_old = V_nxt;
    end
    
    for i = 1:nkgrid
        for j = 1:nkgrid
          rtrn = A*(k_grid(i)^alpha) + (1-delta)*k_grid(i) - kp_grid(j);
          if rtrn < 0
            temp(i, j) = -100.0;
          else
            temp(i, j) = log(rtrn) + beta*V_old(j);
          end
        end
    end

    [argvalue, argmax]= max(temp,[],2);
    dr = argmax;
    V_nxt = argvalue;
    
    for i = 1:nkgrid
        kp_pol_new(i)   = k_grid(argmax(i));
    end
    
    % Apply MacQueen-Porteus Bound in every 5 iteration
    if (mod(iter,5)==0 || iter==1)
        c_low = (beta/(1-beta))*min(V_nxt-V_old);
        c_high = (beta/(1-beta))*max(V_nxt-V_old);
        maxDiff = c_high-c_low;
        fprintf('Current Iteration is: %.d \n', maxDiff);
    end
    
    % Apply MacQueen-Porteus Bound in every iteration
    %c_low = (beta/(1-beta))*min(V_nxt-V_old);
    %c_high = (beta/(1-beta))*max(V_nxt-V_old);
    %maxDiff = c_high-c_low;
    
    V_old = V_nxt;  
    kp_pol = kp_pol_new;
    
    if iter > 100
        break;
    end
end
% Take the median at the end
V_nxt = V_nxt + ((c_high+c_low)/2);
toc; 

figure(1)
plot(k_grid, V_nxt, 'LineWidth',2, 'DisplayName','Value Function');
legend(); grid on;
xlim([min(k_grid) max(k_grid)])
title('Value Function when m = 20')
xlabel('level of Capital Today')
ylabel('Value')

figure(2)
plot(k_grid, kp_pol_new, 'LineWidth',2, 'DisplayName','Capital Policy');
legend(); grid on;
xlim([min(k_grid) max(k_grid)])
title('Optimal Capital Policy Function m = 20')
xlabel('level of Capital Today')
ylabel('level of Capital Tomorrow')
