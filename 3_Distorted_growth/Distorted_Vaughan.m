clear all
clc
tic;

% Parameterization
theta = 0.35; delta = 0.0127; gamma_z = 0.0047; gamma_n = ﻿0.0025;
beta = 0.9937; beta_hat = beta*(1+gamma_n); psi = 2.5;

P =[[ 0.9887,  0.0307, -0.0089,  -0.0407];
    [-0.0012,  1.0011, -0.0275,   0.0175];
    [-0.0045,  0.0449,  0.9675,  -0.0426];
    [ 0.0063,  0.0017,  0.0016,   0.9945]];

P0 = [0.014; 0.0008; 0.0129; -0.0137];

S_bar = inv(eye(4)-P)*P0;
z_ss = exp(S_bar(1));
tl_ss = exp(S_bar(2));
tx_ss = exp(S_bar(3));
g_ss = exp(S_bar(4));

param = [theta, delta, gamma_z, gamma_n, beta_hat, psi];
grate = (1+param(3))*(1+param(4));

%% STEP 1. % Find Steady State
% Find a Steady State
% Kss, Css, Hss obatained analytically.
kl_rto = ( ((1+tx_ss)*(grate-beta_hat*(1-delta))) / (beta_hat*theta*(z_ss^(1-theta))))^(1/(theta-1));
zt1 = (((kl_rto)^(theta-1))*(z_ss^(1-theta))-grate+1-delta);
zt2 = (1/psi)*((1-tl_ss)*(1-theta)*((kl_rto)^theta)*(z_ss^(1-theta)));
zt3 = zt2/kl_rto;

kss = (zt2+g_ss)/(zt1+zt3);
css = zt1*kss - g_ss;
hss = (1/kl_rto)*kss;
xss = grate*kss - (1-delta)*kss;


tl_ss = S_bar(2);
tx_ss = S_bar(3);
%% STEP 2. % Log-linearize the Intratemporal equation
syms h z tl g kp k 
c = (k^theta)*((z*h)^(1-theta)) - grate*kp + (1-delta)*k - g;
w = (1-theta)*((k/h)^theta)*(z^(1-theta));

l_int = psi*(((k^theta)*((z*h)^(1-theta))) - grate*kp + (1-delta)*k - g) - ((1-h)*(1-exp(tl))*(1-theta)*((k/h)^theta)*(z^(1-theta)));
%l_int = (1-tl)*w*(1-h)-psi*c;
grad_int = gradient(l_int, [k, kp, h, z, tl, g]);
Ix = double(subs(grad_int, [k, kp, h, z, tl, g], [kss, kss, hss, z_ss, tl_ss, g_ss]));

coeff = [Ix(1)*kss;
         Ix(2)*kss;
         Ix(3)*hss;
         Ix(4)*z_ss;
         (Ix(5)/tl_ss)*exp(tl_ss);
         Ix(6)*g_ss];

a1 = coeff(1);
a2 = coeff(2);
a3 = coeff(3);
a4 = coeff(4);
a5 = coeff(5);
a6 = coeff(6);

clear c h z tl g kp k 

%% STEP 3. Log-Linearize Euler Equation
syms k kp kpp h hp z tx g zp txp gp

c = (k^theta)*((z*h)^(1-theta)) - grate*kp + (1-delta)*k - g;
cp = (kp^theta)*((zp*hp)^(1-theta)) - grate*kpp + (1-delta)*kp - gp;
rp = theta*(kp^(theta-1))*((zp*hp)^(1-theta));

l_Eul = -beta_hat*c*(rp + (1-delta)*(1+exp(txp))) + cp*(1+exp(tx))*grate;


grad_Eul = gradient(l_Eul, [k, kp, kpp, h, hp, z, tx, g, zp, txp, gp]);
Ex = double(subs(grad_Eul, [k, kp, kpp, h, hp, z, tx, g, zp, txp, gp], [kss, kss, kss, hss, hss, z_ss, tx_ss, g_ss, z_ss, tx_ss, g_ss]));


% 11 arguments
coeff_Eul = [Ex(1)*kss;
             Ex(2)*kss;
             Ex(3)*kss;
             Ex(4)*hss;
             Ex(5)*hss;
             Ex(6)*z_ss;
             (Ex(7)/tx_ss)*exp(tx_ss);
             Ex(8)*g_ss;
             Ex(9)*z_ss;
             (Ex(10)/tx_ss)*exp(tx_ss);
             Ex(11)*g_ss];

b1 = coeff_Eul(1);
b2 = coeff_Eul(2);
b3 = coeff_Eul(3);
b4 = coeff_Eul(4);
b5 = coeff_Eul(5);
b6 = coeff_Eul(6);
b7 = coeff_Eul(7);
b8 = coeff_Eul(8);
b9 = coeff_Eul(9);
b10 = coeff_Eul(10);
b11 = coeff_Eul(11);

clear k kp kpp h hp z tx g zp txp gp


%% Solve the Generalized Eigenvalue Problem.

A1 = [[1, 0, 0];
      [0, 0, 0];
      [0, b3, b5]];

A2 = [[0, -1, 0];
      [a1, a2, a3];
      [b1, b2, b4]];
  
[eigVec, eigVal] = eig(A2, -A1);

%msg = 'STOP: YOU NEED TO CHECK YOUR ORDER OF EIGENVALUES';
%error(msg)

% Re-ordering if you need.
eigVal = [[eigVal(2,2), 0, 0];
          [0, eigVal(1,1), 0];
          [0,0, eigVal(3,3)]];
      
eigVec = [eigVec(:,2), eigVec(:,1), eigVec(:,3)];

V11 = eigVec(1,1);
inn = eigVal(1,1);
V21 = [eigVec(2,1); eigVec(3,1)];

A = V11*inn*inv(V11);
C = V21*inv(V11);

C1 = C(1);
C2 = C(2);

%% Solve System of Equation.
syms x1 x2 x3 x4 d1 d2 d3 d4
B = [x1, x2, x3, x4];
D2 = [d1, d2, d3, d4];

% 1st Eq
eq1 = a2*B + a3*D2 + [a4, a5, 0, a6];
% 2nd Eq
eq2 = b2*B + b3*A*B + b3*B*P + b4*D2 + b5*C2*B + b5*B*P + [b6, 0, b7, b8] + [b9, 0, b10, b11]*P;

sol = vpasolve([eq1(1)==0, eq1(2)==0, eq1(3)==0, eq1(4)==0,...
                eq2(1)==0, eq2(2)==0, eq2(3)==0, eq2(4)==0], ...
                [x1 x2 x3 x4 d1 d2 d3 d4]);
   
x1_sol = sol.x1;
x2_sol = sol.x2;
x3_sol = sol.x3;
x4_sol = sol.x4;

d1_sol = sol.d1;
d2_sol = sol.d2;
d3_sol = sol.d3;
d4_sol = sol.d4;

B = [x1_sol, x2_sol, x3_sol, x4_sol];
D2 = [d1_sol, d2_sol, d3_sol, d4_sol];

B = double(B);
C = double(C);
D2 = double(D2);
%% Construct policy functions
% Construct Capital Policy function  
k_grid = linspace(0.75*kss, 1.25*kss, 100);
k_pol = zeros(100, 1);

for i = 1:100
    k_pol(i) = A*(log(k_grid(i)/kss)) + B*[0; 1; 1; 0];
end

% Construct Labour Policy function
h_pol = zeros(100, 1);

for i = 1:100
    h_pol(i) = C2*log(k_grid(i)/kss) + D2*[0; 1; 1; 0];
end


figure(1);
tiledlayout(1,2)
nexttile
plot(k_grid, exp(k_pol + log(kss)), 'LineWidth', 2, 'DisplayName','Capital Policy function');
title('Capital Policy Function (Variant of Vaughan)');
yline(kss); xline(kss)
xlabel('level of Capital Today')
ylabel('level of Capital Tomorrow')

nexttile
plot(k_grid, exp(h_pol + log(hss)), 'LineWidth', 2, 'DisplayName','Labour Policy function');
title('Labour Policy Function (Variant of Vaughan)');
xlabel('level of Capital Today')
yline(hss); xline(kss)
ylabel('level of Hours')
saveas(figure(1), 'opt_cap_VV.jpg')


%% Simulate Time Series
% S_t = P*S_t-1 + Q*eps

Q =[[ 0.0077,  0.0024, -0.0041,  0.0003];
    [ 0.0024,  0.0043,  0.0023,  0.0153];
    [-0.0041,  0.0023,  0.0088,  0.0121];
    [ 0.0003,  0.0153,  0.0121,  0.0139]];

V = Q*Q';
% Simulate for 1000 periods
horiz = 2000;
St = zeros(4,horiz);

% Setup starting point
S0 = mvnrnd([0,0,0,0], V)';
eps0 = mvnrnd([0,0,0,0], V);
eps0 = eps0';
St(:,1) = P*S0 +  eps0;


for i = 2:horiz
    eps = mvnrnd([0 0 0 0], V);
    St(:,i) =  P*St(:,i-1)  + eps';
end

timesp = linspace(1, horiz, horiz);
z_trend = zeros(horiz, 1);
tl_trend = zeros(horiz, 1);
tx_trend = zeros(horiz, 1);
g_trend = zeros(horiz, 1);

for i = 1:horiz
    z_trend(i) = z_ss*exp(St(1,i));
    tl_trend(i) = tl_ss*(St(2,i));
    tx_trend(i) = tx_ss*(St(3,i));
    g_trend(i) = g_ss*exp(St(4,i));
end

% Construct Time series for Wedges first
figure(2)
plot(timesp, z_trend, 'LineWidth',2, 'DisplayName','Productivity')
title('Simulation on Wedges');
xlabel('Periods')
ylabel('level')
hold on
plot(timesp, tl_trend, 'LineWidth',2, 'DisplayName','Labour Tax')
plot(timesp, tx_trend, 'LineWidth',2, 'DisplayName','Invest Tax')
plot(timesp, g_trend, 'LineWidth',2, 'DisplayName','Gov Exp')
legend;
hold off
saveas(figure(2), 'simulated_wedge.jpg')

% Construct Time Series for Capital and Labor
kstate = zeros(horiz, 1);
lstate = zeros(horiz, 1);
k0 = 2;

kstate(1) = exp((A*log(k0/kss)) + B*S0 +  log(kss));
lstate(1) = exp((C2*log(k0/kss)) + D2*S0 + log(hss));

for i = 2:horiz
    kstate(i) = exp(A*log(kstate(i-1)/kss) + B*St(:,i-1) + log(kss)); 
    lstate(i) = exp(C2*log(kstate(i)/kss) + D2*St(:,i) + log(hss));
end


figure(3)
plot(timesp, kstate, 'LineWidth',2, 'DisplayName','Capital');
title('Simulation of Capital and Labour over time');
xlabel('Periods')
ylabel('level')
yline(kss, 'DisplayName','Capital SS')
hold on
plot(timesp, lstate, 'LineWidth',2, 'DisplayName','Hours');
yline(hss, 'DisplayName','Hours SS')
legend()
hold off
saveas(figure(3), 'sim_cap_hrs.jpg')


% Export result of A, B, C2, D2 from for Next step.
writematrix(A,'A.txt','Delimiter','tab')
writematrix(B,'B.txt','Delimiter','tab')
writematrix(C,'C.txt','Delimiter','tab')
writematrix(D2,'D2.txt','Delimiter','tab')


%% Part 2. Construct Dividend, Profit, and Stock Valuation
% x(t) = grate*k(t+1) - (1-delta)*k(t)

% I simply used 50 periods
horiz = 35;
xstate = zeros(1, horiz-1);
ystate = zeros(1, horiz-1);
cstate = zeros(1, horiz-1);
wage = zeros(1, horiz-1);

for i = 1:horiz-1
    xstate(i) = grate*kstate(i+1) - (1-delta)*kstate(i);
    ystate(i) = (kstate(i)^theta)*((z_trend(i)*lstate(i))^(1-theta));
    cstate(i) = ystate(i) - xstate(i) - g_trend(i);
    wage(i) = z_trend(i)*(1-theta)*(kstate(i)^theta)*((z_trend(i)*lstate(i))^-theta);
end
% Construct dividend and accounting profits
divi = zeros(1, horiz-1);
prft = zeros(1, horiz-1);
Ystate = zeros(1, horiz-1);
for i = 1:horiz-1
    divi(i) = ((grate)^i)*(ystate(i) - (1-tl_trend(i))*wage(i)*lstate(i) - (1-tx_trend(i))*xstate(i));
    prft(i) = ((grate)^i)*(ystate(i) - (1-tl_trend(i))*wage(i)*lstate(i) - (1-tx_trend(i))*delta*kstate(i));
    Ystate(i) = ((grate)^i)*ystate(i);
end

% Construct Arrow-Debreu Price
adp = zeros(1, horiz-1);
adp(1) = 1;
for i = 1:horiz-2
    adp(i+1) = beta*(1+gamma_n)*(cstate(i)/(cstate(i+1)*(1+gamma_z)))*adp(i);
end
SDF = zeros(1, horiz-1);
ctemp = (kstate(35)^theta)*((z_trend(35)*lstate(35))^(1-theta)) - (grate*kstate(36) - (1-delta)*kstate(35)) - g_trend(35);
SDF(34) = beta*(1+gamma_n)*(cstate(34)/(ctemp*(1+gamma_z)));
for i = 1:horiz-2
    SDF(i) = beta*(1+gamma_n)*(cstate(i)/(cstate(i+1)*(1+gamma_z)));
end

% Construct Alternative Price
%DJ_price = readmatrix('DJ_price_ind.csv');
%DJ_price = DJ_price(1:50, 1);
%DJ_price = DJ_price';

% Construct Stock Valuation
stck_val = adp.*divi;

stck_relY = stck_val./Ystate;

figure(4)
timesp3 = linspace(1, horiz-1, horiz-1);
plot(timesp3, divi, 'LineWidth',2, 'DisplayName','Total Dividend');
title('Aggregate Total Net Dividend, Accounting Profits, and Stock Valuation');
xlabel('Periods');
ylabel('Level');
hold on;
plot(timesp3, prft, 'LineWidth',2, 'DisplayName','Accounting Profits');
plot(timesp3, stck_val, 'LineWidth',2, 'DisplayName','Stock Valuation');
legend('Location','northwest');
hold off;

figure(5)
timesp3 = linspace(1, horiz-1, horiz-1);
plot(timesp3, SDF, 'LineWidth',2, 'DisplayName','Stochastic Discount Factor');
title('Stochastic Discount Factor of this Economy');
xlabel('Periods');
ylabel('Level');
legend('location','best');

mean(SDF);
std(SDF);

% Construct Alternative Price
ex_rtrn = readmatrix('excess_rtrn.csv');
ex_rtrn = ex_rtrn(:,2);
mean(ex_rtrn);
std(ex_rtrn);

srp_rtio = mean(ex_rtrn)/std(ex_rtrn);

toc;
