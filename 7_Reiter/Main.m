clear all; clc;
tic;

% Initialization
addpath('lib_ha'); addpath('lib_mirandafackler'); addpath('gensys');
global Params; initparams;

%% STEP 1: Construct Grids on K and Distribution
% Define a grid over cross-sectional distribution.
ndk = Params.ndstst-1;
[Params.knotDistrK, Params.logshift] = makeknotd(Params.kmin, Params.kmax, ndk);

% Define Knot point for Saving polynomial
xmin = Params.kmin;
xmax = Params.kmax; 
n1 = 15;
x1 = linspace(xmin,0.2,n1+1)';
x2 = logspaceshift(0.2, xmax, Params.n_k - n1, Params.logshift)';
knotXi = [x1(1:end-1);x2];
knotXi = knotXi(2:end);
Params.k_grid = [0.01; knotXi];


%% STEP 2: Compute steady state
[Xss, Rss, Wss] = find_ststate;

%% STEP 3: Construct linear system matrices
% Input the system form: G0*X(t) = G1*X(t-1) + C + Psi*z(t) + Pi*eta(t)
% Find G0 and G1 by linearing euler residual equation.
Params.n_states = length(Xss);
Fx = @(X, XP) res_equation(X, XP);

[G0, G1] = linearize(Fx, Xss, Params.epsilon);

% Fill in the C, Psi, and Pi.
Psi = zeros(Params.n_states, 1); Psi(2059) = 1;
Pi = [zeros(Params.n_states - 60, 60); eye(60)];
C  =  zeros(Params.n_states, 1);
%% STEP 4: Solve the linear Rational Expectations system using gensys
%  X(t) = A*X(t-1) + Ctilde + B*z(t) + ywt*inv(I-fmat*inv(L))*fwt*z(t+1) .
[A, ~, B, ~, ~, ~, ~, eu] = gensys(G0, G1, C, Psi, Pi);

fprintf('--------------------------------\n')
fprintf('Existence of eqm: \t %i\n', eu(1))
fprintf('Uniqueness of eqm: \t %i\n', eu(2))
fprintf('--------------------------------\n')

writematrix(A, 'A_mat.txt', 'Delimiter', 'tab');
writematrix(B, 'B_mat.txt', 'Delimiter', 'tab');

A = readmatrix('A_mat.txt');
B = readmatrix('B_mat.txt');

toc;
%% Graphes
% Plot steady state policy functions and distribution
kp_star_ss = reshape(Xss(1:60), Params.n_k, Params.n_y);
hist_ss = reshape(Xss(61:2058), Params.ndstst-1, Params.n_y);
hist_ss_pdf = Xss(61:1059) + Xss(1060:2058);
K_ss = Params.knotDistrK'*hist_ss_pdf;

% Steady State Consumption
for i = 1:Params.n_k
    consm(i, 1) = Rss*Params.k_grid(i) - kp_star_ss(i, 1);
    consm(i, 2) = Rss*Params.k_grid(i) + Wss - kp_star_ss(i, 2);
end

figure(1);
subplot(1,2,1);
legend;
plot(Params.k_grid, kp_star_ss, 'linewidth', 1.5);
grid on;
title('Capital Policy function at Steady State');
xlabel('Capital Today'); ylabel('Capital Tomorrow');
legend('Low productivity','High productivity','Location','Best')

subplot(1,2,2)
plot(Params.k_grid, consm, 'linewidth', 1.5)
grid on
title('Consumption Policy function at Steady State')
xlabel('Capital Today')
ylabel('Consumption Today')
legend('Low productivity','High productivity','Location','northwest')

figure(2)
plot(Params.knotDistrK, hist_ss, 'linewidth', 2)
grid on
title('Distribution of capital in SS')
xlabel('Current Capital Level')
ylabel('Density'); xlim([1, 4])
legend('Low productivity','High productivity','Location','Best')


% Simulating using A and B
% Plot impulse response functions to aggregate shock
T = 500;
simul = simult(A, B, Xss, T);
figure(3)
subplot(3,1,1);
plot(linspace(1,T,T), simul(1,:), 'linewidth', 1.5, 'DisplayName', 'Aggregate Capital') 
grid on; title('Equilibrium Aggregate Capital over time'); xlabel('Time');
yline(K_ss, 'DisplayName', 'Kss')
ylabel('Aggregate Capital'); legend();
subplot(3,1,2);
plot(linspace(1,T,T), simul(2,:), 'linewidth', 1.5, 'DisplayName', 'Gross Interest Rate') 
grid on; title('Equilibrium Gross Interest Rate over time'); xlabel('Time');
yline(Rss, 'DisplayName', 'Rss')
ylabel('R'); legend();
subplot(3,1,3);
plot(linspace(1,T,T), simul(3,:), 'linewidth', 1.5, 'DisplayName', 'Wage') 
grid on; title('Equilibrium Wage over time'); xlabel('Time');
yline(Wss, 'DisplayName', 'Wss')
ylabel('W'); legend();

T = 50;
% Impluse Response
simul_irf = gen_IRF(A, B, T);
for i = 1:T
    K_irf(i) = Params.knotDistrK'*simul_irf(61:1059,i) + Params.knotDistrK'*simul_irf(1060:2058,i);
    R_irf(i) = intr((K_ss+K_irf(i)/0.9300), exp(simul_irf(2059,i)));
    W_irf(i) = cwage((K_ss+K_irf(i)/0.9300), exp(simul_irf(2059,i)));
    
end

figure(4)
subplot(3,1,1);
plot(linspace(1,50,50), K_irf+K_ss, 'linewidth', 1.5, 'DisplayName', 'Aggregate K') 
grid on; title('Impulse Response on K'); xlabel('Time');
ylabel('Aggregate Capital'); legend();
subplot(3,1,2);
plot(linspace(1,50,50), R_irf, 'linewidth', 1.5, 'DisplayName', 'Gross Interest Rate') 
grid on; title('Impulse Response on R'); xlabel('Time');
ylabel('R'); legend();
subplot(3,1,3);
plot(linspace(1,50,50), W_irf, 'linewidth', 1.5, 'DisplayName', 'Wage') 
grid on; title('Impulse Response on W'); xlabel('Time');
ylabel('W'); legend();


