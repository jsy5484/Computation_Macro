clear all clc
tic;

% Parameterization
alpha = 0.3333333333; delta = 0.1; beta = 0.975; A = 1.05;
nkgrid = 1001; maxDiff = 10.0; tol = 1.0e-10; iter = 0;

% Find the Steady-states
kss   = ((1 + beta * (delta - 1))/(alpha * beta * A))^(1 / (alpha - 1));
yss   = A * kss^alpha + (1 - delta)*kss;
css   = yss - kss;

% Construct Capital Grid
kmin = 0.5*kss; kmax = 1.5*kss;
k_grid = linspace(kmin, kmax, nkgrid);
kp_grid = linspace(kmin, kmax, nkgrid);

%rtrn = zeros(ngrid_k, ngrid_k);
V_old = zeros(1, nkgrid);
V_new = zeros(1, nkgrid);
k_pol = zeros(1, nkgrid);

% Main loop
while maxDiff > tol
    iter = iter + 1;
    fprintf('Current Iteration is: %.d \n', iter);
    
    temp = zeros(nkgrid, nkgrid);
    
    % Clear negative consumption first
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
    
    % Find Max and argmax
    [argvalue, argmax]= max(temp,[],2);
    V_new = argvalue;
    
   for i = 1:nkgrid
        k_pol(i) = k_grid(argmax(i));
   end
   
    % Apply MacQueen-Porteus Bound in every iteration
   c_low = (beta/(1-beta))*min(V_new-V_old);
   c_high = (beta/(1-beta))*max(V_new-V_old);
   maxDiff = c_high-c_low;
    
   %maxDiff = norm(V_new-V_old);
   V_old = V_new;  
end
toc; 

V_new = V_new + ((c_high+c_low)/2);


figure(1)
plot(k_grid, V_new, 'LineWidth',2, 'DisplayName','Value Function');
legend(); grid on;
xlim([min(k_grid) max(k_grid)])
title('Value Function Converged')
xlabel('level of Capital Today')
ylabel('Value')

figure(2)
plot(k_grid, k_pol, 'LineWidth',2, 'DisplayName','Capital Policy');
legend(); grid on;
xlim([min(k_grid) max(k_grid)])
title('Optimal Capital Policy Function')
xlabel('level of Capital Today')
ylabel('level of Capital Tomorrow')



