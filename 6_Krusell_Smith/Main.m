clear all; clc;
tic;

% Parameterization
alpha = 0.33; beta = 0.96; delta = 0.05; mu = 3;
param = [alpha beta delta mu];
% Unemployment rate 
unemp_gd = 0.04; unemp_bd = 0.1;
% Simulation Horizon
horiz = 1000;
% length of 
dur_ug = 1.5; dur_ub = 2.5; dur_gd = 8.0; dur_bd = 8.0;
% Shocks
eps_gd = 1; eps_bd = 0; z_gd = 1.01; z_bd = 0.99;

%% Step 1. Construct a Transition Matrix 
Trans = constructTrans(dur_ug, dur_ub, dur_gd, dur_bd, unemp_gd, unemp_bd);
Trans = [Trans(4,:);Trans(2,:);Trans(3,:);Trans(1,:)];
Trans = [Trans(:,4),Trans(:,2),Trans(:,3),Trans(:,1)];

% Get simulation about shocks, assumming 10,000 household
N = 10000;
[agg_shock, idio_shock] = shock_simul(Trans, N, horiz, unemp_bd);     

%% Step 2. Construct Grids.
% Grid for the first moment of Aggregate
KM_min = 40; KM_max = 120; ngrid_KM = 100;
KM_grid = linspace(KM_min, KM_max, ngrid_KM)';

% Grid for the individual
k_min = 0; k_max = 400; ngrid_k = 1000;
x = linspace(0, 0.5, ngrid_k)';
y = x.^7/max(x.^7);
k_grid = k_min + (k_max - k_min)*y;
clear x y

% Fill in the initial level of Aggregate Capital and individual Capital
% (1) Initial guess of regression coefficient, B. 
B = [0, 1, 0, 1];

% (2) Initial guess for policy function.
k_pol = zeros(ngrid_k, ngrid_KM, 2, 2);
for i = 1:ngrid_KM
    for j = 1:2
        for w = 1:2
            k_pol(:,i,j,w) = 0.9*k_grid;
        end
    end
end
% (3) Initial guess for Cross-Sectional Distribution : Steady-State Level
kss = ((1/beta-(1-delta))/alpha)^(1/(alpha-1));
KM_cross = zeros(1, N) + kss;

%% STEP 3. Main Loop
iter = 0; diff = 100; tol = 1e-8; B1 = [0, 0, 0, 0];

while diff > tol
    iter = iter + 1;
    fprintf('Current Iteration is : %d \n', iter);
    
    % Solve Individual Problem using Policy function Iteration.
    [k_pol, consm] = solveIND(Trans, unemp_gd, unemp_bd, ngrid_k, ngrid_KM, k_grid, KM_grid, k_pol, param, B, KM_min, KM_max, k_min, k_max);
    
    % Do simulation with obtained optimal policy function, 
    % return time-series simulation and cross-sectional distribution.
    [KM_time, KM_cross_new] = aggSimul(horiz, idio_shock, agg_shock, KM_max, KM_min, k_pol, KM_grid, k_grid, k_min, k_max, KM_cross);
    
    % Regression Part / Set a counter to capture good/bad time.
    bd_counter = 0; gd_counter = 0;
    x_var_bd = 0; y_var_bd = 0; x_var_gd = 0; y_var_gd = 0;
    
    for i = 1:(horiz-1)
        if agg_shock(i) == 1
            bd_counter = bd_counter+1;
            x_var_bd(bd_counter, 1) = log(KM_time(i));
            y_var_bd(bd_counter, 1) = log(KM_time(i+1));
        else
            gd_counter = gd_counter+1;
            x_var_gd(gd_counter, 1) = log(KM_time(i));
            y_var_gd(gd_counter, 1) = log(KM_time(i+1));
        end
    end
    
    [B1(1:2), sb2, sb3, sb4, s5] = regress(y_var_bd, [ones(bd_counter, 1), x_var_bd]);
    r_sqrd_bd = s5(1);
    
    [B1(3:4), sg2, sg3, sg4, s5] = regress(y_var_gd, [ones(gd_counter, 1), x_var_gd]);
    r_sqrd_gd = s5(1);
    
    diff = norm(B - B1);
    fprintf('R_squared difference is : %d \n', diff);
    
    if diff >(tol*100)
        KM_cross = KM_cross_new; 
    end

    B = B1;                            
end

% Construct time-series of Aggregate Capital using obtained B1, B2 to
% compare
KM_sol = zeros(N,1);  
KM_sol(1)= KM_time(1); 
for i = 1:horiz-1       
   if agg_shock(i)==1
      KM_sol(i+1)=exp(B1(1)+B1(2)*log(KM_sol(i)));
   else
      KM_sol(i+1)=exp(B1(3)+B1(4)*log(KM_sol(i)));
   end
end



%% Figure
figure(1)
timesp = linspace(1,1000,1000);
plot(timesp, KM_time(1:horiz,1), 'LineWidth',1, 'DisplayName','Agg Capital implied by individual policy');
title('Comparison of Actual and Estimated Law of Motion')
hold on; grid on
plot(timesp, KM_sol(1:horiz,1), 'LineWidth',1, 'DisplayName','Agg Law of Motion for Captital');
legend;
xlabel('Time'), ylabel('Aggregate capital')
%saveas(figure(1), 'comparision.jpg')

figure(2)
timesp1 = linspace(1, 100, 100);
plot(timesp1, exp(B1(1)+B1(2)*log(KM_grid)), 'LineWidth',1, 'DisplayName','Bad');
title('Law of Motion for Aggregate Capital');
hold on; grid on;
plot(timesp1, exp(B1(3)+B1(4)*log(KM_grid)), 'LineWidth',1, 'DisplayName','Good');
legend; xlim([20 30]); 
xlabel('Aggregate Capital Today'),  ylabel('Aggregate Capital Tomorrow')
%saveas(figure(2), 'AggLM.jpg')


toc;


