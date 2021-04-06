clear all
clc; tic;

global Params
initparams; 
tol = 1.0e-6; iter = 0; maxDiff = 100;

% Find the Steady-state
kss   = ((1 + Params.beta * (Params.delta - 1))/(Params.theta * Params.beta))^(1 / (Params.theta - 1));
yss   = kss^Params.theta + (1 - Params.delta)*kss;
css   = yss - kss;

% Discretize the stochastic process
[Params.trans, z_grid] = rouwTrans(Params.rho, 0, Params.sigma/sqrt(Params.ngrid_z-1), Params.ngrid_z);
%[Params.trans, z_grid] = rouwen(Params.rho, 0, Params.sigma/sqrt(1-Params.rho^2), Params.ngrid_z);
Params.z_grid = exp(z_grid);
% Construct Capital Grid
Params.kmin = z_grid(1)*(-0.6);
Params.kmax = 10;
%Params.k_grid = linspace(Params.kmin, Params.kmax, Params.ngrid_k);
Params.k_grid = Params.kmin + (Params.kmax - Params.kmin)*linspace(0, 1, Params.ngrid_k).^3;


k_grid = Params.k_grid';
kp_grid = Params.k_grid';

C_next  = zeros(Params.ngrid_k, Params.ngrid_z);
C_star  = zeros(Params.ngrid_k, Params.ngrid_z);
B_pol   = zeros(Params.ngrid_k, Params.ngrid_z);

% A guess for consumption
for j = 1:Params.ngrid_z
    for i = 1:Params.ngrid_k
        Y_end(i, j) = (Params.z_grid(j)*(Params.k_grid(i)^Params.theta) + (1 - Params.delta)*Params.k_grid(i));
        util = (Y_end(i, j)^(1-Params.alpha))/(1-Params.alpha);
        V_old(i, j) = util/(1-Params.beta);
    end
end


while (maxDiff > tol)
    fprintf('Current iteration is: %.d \n', iter);
    fprintf('Current Difference is: %.8f \n', maxDiff);

    iter = iter + 1;
    
    % Find the Envelope condition using Numerical difference
    V_env = zeros(Params.ngrid_k, Params.ngrid_z); 
    V_env(1,:) = ((V_old(2,:)-V_old(1,:))/(k_grid(2)-k_grid(1)));
    V_env(Params.ngrid_k,:) = ((V_old(Params.ngrid_k,:)-V_old(Params.ngrid_k-1,:))/(k_grid(Params.ngrid_k)-k_grid(Params.ngrid_k-1)));
    for j = 1:Params.ngrid_z
        for i = 2:Params.ngrid_k-1
            V_env(i,j) = (0.5)*((V_old(i,j)-V_old(i-1,j))/(k_grid(i)-k_grid(i-1))) + (0.5)*((V_old(i+1,j)-V_old(i,j))/(k_grid(i+1)-k_grid(i)));
        end
    end
    
    
    for j = 1:Params.ngrid_z
        for i = 1:Params.ngrid_k
            C_star(i, j) = V_env(i,j)^(-1/Params.alpha);
            Y_star(i, j) = C_star(i, j) + Params.k_grid(i);
            V_int(i, j) = ((C_star(i,j)^(1-Params.alpha))/(1-Params.alpha)) + V_old(i,j);
        end
    end
    
    % Update the Value Function
    %V_int = ((C_star.^(1-Params.alpha))/(1-Params.alpha)) + V_old;
    for j = 1:Params.ngrid_z
        for i = 1:Params.ngrid_k
            V_end(i,j) = interp1(Y_star(:,j), V_int(:,j), Y_end(i,j), 'spline', 'extrap');
            %V_end(i,j) = interp1(Y_end(:,j), V_int(:,j), Y_star(i,j), 'linear', 'extrap');
        end
    end
    
    V_nxt = Params.beta*V_end*Params.trans;
    
    maxDiff = norm(V_nxt - V_old);
    V_old = V_nxt;
    
end

for j = 1:Params.ngrid_z
    for i = 1:Params.ngrid_k  
        fun = @(x) Params.z_grid(j)*(x^Params.theta) + (1-Params.delta)*x - Y_star(i, j);
        k_back(i, j) = fsolve(fun, 0.1);
        %k_back(i, j) = bisection(fun, -0.1, 10);
    end
end
toc

figure(1)
plot(k_grid, V_nxt(:, 1), 'LineWidth',1.5, 'DisplayName','Low');
hold on;
plot(k_grid, V_nxt(:, 6), 'LineWidth',1.5, 'DisplayName','Mid');
plot(k_grid, V_nxt(:, 11), 'LineWidth',1.5, 'DisplayName','High');
legend('Location','best'); grid on; 
title('Value Function for Low, Mid, High productivity'); xlabel('level of Capital Today');
ylabel('Value'); hold off;

figure(2)
plot(k_back(:, 1), k_grid, 'LineWidth',1.5, 'DisplayName','Low');
hold on; 
plot(k_back(:, 6), k_grid, 'LineWidth',1.5, 'DisplayName','Mid');
plot(k_back(:, 11), k_grid, 'LineWidth',1.5, 'DisplayName','High');
legend('Location','best');  grid on; 
title('Optimal Capital Policy Function for Low, Mid, and High productivity')
xlabel('level of Capital Today'); ylabel('level of Capital Tomorrow')
hold off;