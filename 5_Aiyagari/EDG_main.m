clear all;
clc; tic;

% Parameterization
alpha = 1/3; delta = 0.1; beta = 0.96; mu = 2; l_low = 0.2; l_high = 1.2;
param = [beta, mu, l_low, l_high];

%% Main loop: Find Equilibrium price
tol = 1e-5; iter = 0; diff = 100;

% Initial guess
r_old = 0.03; w_old = (1-alpha)*((r_old+delta)/alpha)^(alpha/(alpha-1));
while  diff > tol
    iter = iter + 1;
    fprintf('Iterating on prices, current iteration is: %.d \n', iter);

    [A, W, state] = EDG_grid(r_old, w_old, param);
    K = sum(W)/10000;
    L = (sum(state)/10000)*l_high + ((1-sum(state))/10000)*l_low;
    
    % demand side
    r_new = alpha*((K/L)^(alpha-1)) - delta;
    w_new = (1-alpha)*((r_new+delta)/alpha)^(alpha/(alpha-1));
    
    
    diff = abs(r_new - r_old);
    fprintf('Interest rate difference is: %.5f \n', diff);
    
    r_old = r_new;
end

R_eq = r_new+1; 
w_eq = w_new;

toc;

fprintf('Equilibrium Interest Rate is : %.6f \n', R_eq);
fprintf('Equilibrium Wage is: %.4f \n', w_eq);
