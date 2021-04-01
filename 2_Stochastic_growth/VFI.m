clear all
clc
tic;

theta = 0.35; delta = 0.0464; gamma_z = 0.016; gamma_n = 0.015;
beta = 0.9722; beta_hat = beta*(1+gamma_n); sigma = 0.5;
rho = 0.2; psi = 2.24;
grate = (1+gamma_z)*(1+gamma_n);
param = [theta, delta, gamma_z, gamma_n, beta_hat, psi];

ngrid_z = 2;
ngrid_k = 1000;

% Get the Steady-State level of capital from a Nonlinear solver.
SS = SS(param);
kss = SS(1);
hss = SS(2);
lss = 1 - hss;
css = ((kss^param(1))*((exp(0)*hss)^(1-param(1)))) - (1+param(4))*(1+param(3))*kss + (1-param(2))*kss;

% Discretize the stochastic process
[z_grid, trans] = tauchenHussey(ngrid_z, 0, rho, sigma, sigma);

% Construct Capital Grid
kmin = 0.5*kss;
kmax = 1.5*kss;
k_grid = linspace(kmin, kmax, ngrid_k);
k_grid = k_grid';
kt_grid = k_grid';


rtrn = zeros(ngrid_k, ngrid_k, ngrid_z);
h_star = zeros(ngrid_k, ngrid_k, ngrid_z);

% Solve h from intratemporal equation over the grids
for i = 1:ngrid_k
    for j = 1:ngrid_k
        for w = 1:ngrid_z
            fun = @(hr) (psi/(1-theta))*((k_grid(i)^theta)*(exp(z_grid(w))*hr)^(1-theta) - grate*k_grid(j) + (1-delta)*k_grid(i)) + ((hr-1)*(k_grid(i)^theta)*(exp(z_grid(w))^(1-theta))*(hr^(-theta)));
            %h = fsolve(fun, 0.5);
            h = bisection(fun, -10, 10);
            if h > 1
                h_star(i,j,w) = 0.9999;
            else
                h_star(i,j,w) = h;
            end
        end
    end
end

% Precalculate return value
for i = 1:ngrid_k
    for j = 1:ngrid_k
        for w = 1:ngrid_z
           pre_cnsm = (k_grid(i)^theta)*((exp(z_grid(w))*h_star(i,j,w))^(1-theta))-grate*k_grid(j)+(1-delta)*k_grid(i);
           pre_rttn = log(pre_cnsm) + psi*log(1-h_star(i,j,w));
           
           if pre_cnsm < 0
                rtrn(i,j,w) = -100.0;
           else
                rtrn(i,j,w) = pre_rttn;
           end
       
        end
    end
end

V_old = zeros(ngrid_k, ngrid_z);
% Fill an initial guess
for j = 1:ngrid_z
    for i = 1:ngrid_k
        V_old(i, j) = (log((k_grid(i)^theta)*((exp(z_grid(j))*hss)^(1-theta))-delta*k_grid(i))+psi*log(1-hss))/(1-beta_hat);
    end
end

V_new = zeros(ngrid_k, ngrid_z);
opt_ind = zeros(ngrid_k, ngrid_z);
expectedVal = zeros(ngrid_k, ngrid_z);

maxDiff = 10.0; tol = 1.0e-10; iter = 0;

% Main loop
while maxDiff > tol
    iter = iter + 1;
    fprintf('Current Iteration is: %.d \n', iter);
    
    k_pol = zeros(ngrid_k, ngrid_z);
    h_pol = zeros(ngrid_k, ngrid_z);
    c_pol = zeros(ngrid_k, ngrid_z);

    expected_val = V_old*trans';
    temp = zeros(size(rtrn));
    for j = 1:size(rtrn,1)
        temp(j,:,:) = rtrn(j,:,:) + reshape(beta_hat*expected_val,1,size(expected_val,1),size(expected_val,2));
    end
    [argvalue, argmax]= max(temp,[],2);
    V_new = squeeze(argvalue);
    
    
   for i = 1:ngrid_k
        for j = 1:ngrid_z
          %  [V_new(i,j), opt_ind(i,j)] = max(squeeze(rtrn(i,:,j)) + beta_hat*expected_val(i,j)); 
            k_pol(i,j) = k_grid(argmax(i,j));
            h_pol(i,j) = h_star(i, argmax(i,j), j);
            c_pol(i,j) = (k_grid(i)^theta)*((exp(z_grid(j))*h_pol(i,j))^(1-theta))-grate*k_pol(i,j)+(1-delta)*k_grid(i);
        end
   end
    maxDiff = norm(V_new-V_old);
    V_old = V_new;  
end

figure(1)
plot(k_grid, V_new(:, 1), 'LineWidth',2, 'DisplayName','Low');
hold on;
plot(k_grid, V_new(:, 2), 'LineWidth',2, 'DisplayName','High');
legend();
xlim([kmin kmax])
title('Value Function')
xlabel('level of Capital Today')
ylabel('Value')
hold off

figure(2)
plot(k_grid, k_pol(:, 1), 'LineWidth',2, 'DisplayName','Low');
hold on;
plot(k_grid, k_pol(:, 2), 'LineWidth',2, 'DisplayName','High');
legend();
xlim([kmin kmax])
title('Optimal Capital Policy Function')
xlabel('level of Capital Today')
ylabel('level of Capital Tomorrow')
hold off

figure(3)
plot(k_grid, h_pol(:, 1), 'LineWidth',2, 'DisplayName','Low');
hold on;
plot(k_grid, h_pol(:, 2), 'LineWidth',2, 'DisplayName','High');
xlim([kmin kmax])
legend();
title('Optimal Labour Policy Function')
xlabel('level of Capital Today')
ylabel('level of Hours')
hold off

figure(4)
plot(k_grid, c_pol(:, 1), 'LineWidth',2, 'DisplayName','Low');
hold on;
plot(k_grid, c_pol(:, 2), 'LineWidth',2, 'DisplayName','High');
legend();
xlim([kmin kmax])
title('Optimal Consumption')
xlabel('level of Capital Today')
ylabel('level of Consumption')
hold off

toc;