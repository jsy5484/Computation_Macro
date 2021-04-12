clear all
clc; tic;

global Params
initparams; 

% Find the Steady-state
kss   = ((1 + Params.beta * (Params.delta - 1))/(Params.theta * Params.beta))^(1 / (Params.theta - 1));
yss   = kss^Params.theta + (1 - Params.delta)*kss;
css   = yss - kss;

% Discretize the stochastic process
[Params.trans, z_grid] = rouwTrans(Params.rho, 0, Params.sigma/sqrt(Params.ngrid_z-1), Params.ngrid_z);
Params.z_grid = exp(z_grid);

% Construct Capital Grid
Params.kmin = 0.5*kss;
Params.kmax = 1.5*kss;
Params.k_grid = linspace(Params.kmin, Params.kmax, Params.ngrid_k);
k_grid = Params.k_grid';
kp_grid = Params.k_grid';

rtrn = zeros(Params.ngrid_k, Params.ngrid_k, Params.ngrid_z);

% Precalculate return value
for i = 1:Params.ngrid_k
    for j = 1:Params.ngrid_k
        for w = 1:Params.ngrid_z
           pre_cnsm = Params.z_grid(w)*(Params.k_grid(i)^Params.theta) + (1 - Params.delta)*Params.k_grid(i) - kp_grid(j);
           if pre_cnsm < 0
                rtrn(i,j,w) = -3000.0;
           else
                rtrn(i,j,w) = (pre_cnsm^(1-Params.alpha))/(1-Params.alpha);
           end
        end
    end
end

Params.V_old = zeros(Params.ngrid_k, Params.ngrid_z);
% Fill an initial guess
for i = 1:Params.ngrid_k
    for j = 1:Params.ngrid_z
        consm = (Params.z_grid(j)*(Params.k_grid(i)^Params.theta) - Params.delta*Params.k_grid(i));
        util = (consm^(1-Params.alpha))/(1-Params.alpha);
        Params.V_old(i, j) = util/(1-Params.beta);
    end
end

Params.V_new = zeros(Params.ngrid_k, Params.ngrid_z);
expectedVal = zeros(Params.ngrid_k, Params.ngrid_z);

maxDiff = 10.0; tol = 1.0e-10; iter = 0;
while maxDiff > tol
    iter = iter + 1;
    fprintf('Current Iteration is: %.d \n', iter);
    fprintf('Current Iteration is: %.8f \n', maxDiff);

    for j = 1:Params.ngrid_z
        for i = 1:Params.ngrid_k   
            Params.k_now = Params.k_grid(i);
            Params.z_now = Params.z_grid(j);
            Params.temp = j;
            k_mnmizer = fminbnd( @Val_func, Params.kmin, Params.kmax); 
            Params.V_new(i,j) = - Val_func(k_mnmizer);
            k_pol(i,j) = k_mnmizer;
        end
    end
    maxDiff = norm(Params.V_new - Params.V_old);
    Params.V_old = Params.V_new;
end
toc;

figure(1)
plot(k_grid, Params.V_new(:, 1), 'LineWidth',1.5, 'DisplayName','Low');
hold on;
plot(k_grid, Params.V_new(:, 6), 'LineWidth',1.5, 'DisplayName','Mid');
plot(k_grid, Params.V_new(:, 11), 'LineWidth',1.5, 'DisplayName','High');
legend('Location','best'); grid on; xlim([Params.kmin Params.kmax])
title('Value Function for Low, Mid, High productivity'); xlabel('level of Capital Today');
ylabel('Value')
hold off

figure(2)
plot(k_grid, k_pol(:, 1), 'LineWidth',1.5, 'DisplayName','Low');
hold on; 
plot(k_grid, k_pol(:, 6), 'LineWidth',1.5, 'DisplayName','Low');
plot(k_grid, k_pol(:, 11), 'LineWidth',1.5, 'DisplayName','High');
legend('Location','best'); xlim([3, 5]); grid on;
title('Optimal Capital Policy Function for Low, Mid, and High productivity')
xlabel('level of Capital Today'); ylabel('level of Capital Tomorrow')
hold off
