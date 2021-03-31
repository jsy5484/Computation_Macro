% I chose the same parameter as before
clear all;
clc; tic;

theta = 0.35; delta = 0.0464; gamma_z = 0.016; gamma_n = 0.015;
beta = 0.9722; beta_hat = beta*(1+gamma_n); psi = 2.24;
grate = (1+gamma_z)*(1+gamma_n);
P0 = [0.014; 0.0008; 0.0129; -0.0137];

P =[[ 0.9887,  0.0307, -0.0089,  -0.0407];
    [-0.0012,  1.0011, -0.0275,   0.0175];
    [-0.0045,  0.0449,  0.9675,  -0.0426];
    [ 0.0063,  0.0017,  0.0016,   0.9945]];

S_bar = (eye(4)-P)\P0;

Q =[[ 0.0077,  0.0024, -0.0041,  0.0003];
    [ 0.0024,  0.0043,  0.0023,  0.0153];
    [-0.0041,  0.0023,  0.0088,  0.0121];
    [ 0.0003,  0.0153,  0.0121,  0.0139]];
V = Q*Q';

% Steady-States are
z_ss = exp(S_bar(1));
tl_ss = exp(S_bar(2));
tx_ss = exp(S_bar(3));
g_ss = exp(S_bar(4));

kl_rto = ( ((1+tx_ss)*(grate-beta_hat*(1-delta))) / (beta_hat*theta*(z_ss^(1-theta))))^(1/(theta-1));
zt1 = (((kl_rto)^(theta-1))*(z_ss^(1-theta))-grate+1-delta);
zt2 = (1/psi)*((1-tl_ss)*(1-theta)*((kl_rto)^theta)*(z_ss^(1-theta)));
zt3 = zt2/kl_rto;

kss = (zt2+g_ss)/(zt1+zt3);
css = zt1*kss - g_ss;
hss = (1/kl_rto)*kss;
xss = grate*kss - (1-delta)*kss;
yss = (kss^theta)*((z_ss*hss)^(1-theta));

% Construct policy function. Import key ingredient from previous Vaughan's Method
A = readmatrix('A.txt');
B = importdata('B.txt');
C = importdata('C.txt');
C2 = C(2);
D2 = importdata('D2.txt');

%% PART A. Simulation
% N = 200
% Simulate for 200 periods
horiz = 200;
St = zeros(4,horiz);

% Setup starting point
S0 = S_bar';
eps = mvnrnd([0,0,0,0], eye(4));

St(:,1) = P*S0' + Q*eps';

for i = 2:horiz
    eps = mvnrnd([0 0 0 0], eye(4));
    St(:,i) = P*St(:,i-1) + Q*eps';
end

timesp = linspace(1, 200, 200);
z_trend = zeros(200, 1);
tl_trend = zeros(200, 1);
tx_trend = zeros(200, 1);
g_trend = zeros(200, 1);

for i = 1:horiz
    z_trend(i) = (St(1,i));
    tl_trend(i) = (St(2,i));
    tx_trend(i) = (St(3,i));
    g_trend(i) = (St(4,i));
end

kstate = zeros(201, 1);
lstate = zeros(200, 1);
cstate = zeros(200, 1);
ystate = zeros(200, 1);
istate = zeros(200, 1);
k0 = 1.8;

kstate(1) = exp((A*log(k0/kss) + B*S0') + log(kss));

%exp(z_trend(i)-log(z_ss) +   )
for i = 2:201
    kstate(i) = exp((A*log(kstate(i-1)/kss) + B*St(:,i-1)) + log(kss)); 
end

for i = 1:horiz
    lstate(i) = exp((C2*log(kstate(i)/kss) + D2*St(:,i)) + log(hss));
    cstate(i) = (kstate(i)^theta)*((exp(z_trend(i) + log(z_ss))*lstate(i))^(1-theta)) - grate*kstate(i+1) + (1-delta)*kstate(i) - exp(g_trend(i)+log(g_ss));
    istate(i) = grate*kstate(i+1) - (1-delta)*kstate(i);
    ystate(i) = (kstate(i)^theta)*((exp(z_trend(i) + log(z_ss))*lstate(i))^(1-theta)); 
end
kstate(end)=[];

figure(1)
plot(timesp, kstate, 'LineWidth',2, 'DisplayName','Capital')
title('Simulation of C, I, Y, and H for N = 200')
xlabel('Time')
ylabel('level')
hold on
plot(timesp, lstate, 'LineWidth',2, 'DisplayName','Hours');
plot(timesp, istate, 'LineWidth',2, 'DisplayName','Invest');
plot(timesp, cstate, 'LineWidth',2, 'DisplayName','Consumption');
plot(timesp, ystate, 'LineWidth',2, 'DisplayName','Production');
yline(hss, 'DisplayName','Hours SS')
yline(xss, 'DisplayName','Invest SS')
yline(css, 'DisplayName','Consumption SS')
yline(yss, 'DisplayName','Production SS')
yline(kss, 'DisplayName','Captial SS')
legend;
hold off
saveas(figure(1),'simul_200.jpg')

clear z_trend tl_trend tx_trend g_trend kstate lstate istate cstate ystate 

%% N = 1000
% Simulate for 1000 periods
horiz = 1000;
St = zeros(4,horiz);

% Setup starting point
S0 = S_bar;
eps = mvnrnd([0,0,0,0], eye(4));

St(:,1) = P*S0 + Q*eps';

for i = 2:horiz
    eps = mvnrnd([0 0 0 0], eye(4));
    St(:,i) = P*St(:,i-1) + Q*eps';
end


timesp = linspace(1, 1000, 1000);

z_trend = zeros(1000, 1);
tl_trend = zeros(1000, 1);
tx_trend = zeros(1000, 1);
g_trend = zeros(1000, 1);

for i = 1:horiz
    z_trend(i) = (St(1,i));
    tl_trend(i) = (St(2,i));
    tx_trend(i) = (St(3,i));
    g_trend(i) = (St(4,i));
end

kstate = zeros(1001, 1);
lstate = zeros(1000, 1);
cstate = zeros(1000, 1);
ystate = zeros(1000, 1);
istate = zeros(1000, 1);
k0 = 1.8;

kstate(1) = exp((A*log(k0/kss) + B*S0) + log(kss));

%exp(z_trend(i)-log(z_ss) +   )
for i = 2:1001
    kstate(i) = exp((A*log(kstate(i-1)/kss) + B*St(:,i-1)) + log(kss)); 
end

for i = 1:horiz
    lstate(i) = exp((C2*log(kstate(i)/kss) + D2*St(:,i)) + log(hss));
    cstate(i) = (kstate(i)^theta)*((exp(z_trend(i) + log(z_ss))*lstate(i))^(1-theta)) - grate*kstate(i+1) + (1-delta)*kstate(i) - exp(g_trend(i)+log(g_ss));
    istate(i) = grate*kstate(i+1) - (1-delta)*kstate(i);
    ystate(i) = (kstate(i)^theta)*((exp(z_trend(i) + log(z_ss))*lstate(i))^(1-theta)); 
end
kstate(end)=[];

figure(2)
plot(timesp, kstate, 'LineWidth',2, 'DisplayName','Capital')
title('Simulation of C, I, Y, and H for N = 1000')
xlabel('Time')
ylabel('Level')
hold on
plot(timesp, lstate, 'LineWidth',2, 'DisplayName','Hours');
plot(timesp, istate, 'LineWidth',2, 'DisplayName','Invest');
plot(timesp, cstate, 'LineWidth',2, 'DisplayName','Consumption');
plot(timesp, ystate, 'LineWidth',2, 'DisplayName','Production');
yline(hss, 'DisplayName','Hours SS')
yline(xss, 'DisplayName','Invest SS')
yline(css, 'DisplayName','Consumption SS')
yline(yss, 'DisplayName','Production SS')
yline(kss, 'DisplayName','Captial SS')
legend;
hold off
saveas(figure(2),'simul_1000.jpg')

clear z_trend tl_trend tx_trend g_trend kstate lstate istate cstate ystate 

%% N = 5000
% Simulate for 5000 periods
horiz = 5000;
St = zeros(4,horiz);

% Setup starting point
%S0 = mvnrnd([0,0,0,0], eye(4));
S0 = S_bar';
eps = mvnrnd([0,0,0,0], eye(4));

St(:,1) = P*S0' + Q*eps';

for i = 2:horiz
    eps = mvnrnd([0 0 0 0], eye(4));
    St(:,i) = P*St(:,i-1) + Q*eps';
end

timesp = linspace(1, 5000, 5000);

z_trend = zeros(5000, 1);
tl_trend = zeros(5000, 1);
tx_trend = zeros(5000, 1);
g_trend = zeros(5000, 1);

for i = 1:horiz
    z_trend(i) = St(1,i);
    tl_trend(i) = St(2,i);
    tx_trend(i) = St(3,i);
    g_trend(i) = St(4,i);
end

kstate = zeros(5001, 1);
lstate = zeros(5000, 1);
cstate = zeros(5000, 1);
ystate = zeros(5000, 1);
istate = zeros(5000, 1);
k0 = 1.8;

kstate(1) = exp((A*log(k0/kss) + B*S0') + log(kss));

%exp(z_trend(i)-log(z_ss) +   )
for i = 2:5001
    kstate(i) = exp((A*log(kstate(i-1)/kss) + B*St(:,i-1)) + log(kss)); 
end

for i = 1:horiz
    lstate(i) = exp((C2*log(kstate(i)/kss) + D2*St(:,i)) + log(hss));
    cstate(i) = (kstate(i)^theta)*((exp(z_trend(i) + log(z_ss))*lstate(i))^(1-theta)) - grate*kstate(i+1) + (1-delta)*kstate(i) - exp(g_trend(i)+log(g_ss));
    istate(i) = grate*kstate(i+1) - (1-delta)*kstate(i);
    ystate(i) = (kstate(i)^theta)*((exp(z_trend(i) + log(z_ss))*lstate(i))^(1-theta)); 
end
kstate(end)=[];

figure(3)
plot(timesp, kstate, 'LineWidth',2, 'DisplayName','Capital')
title('Simulation of C, I, Y, and H for N = 5000')
xlabel('Time')
ylabel('Level')
hold on
plot(timesp, lstate, 'LineWidth',2, 'DisplayName','Hours');
plot(timesp, istate, 'LineWidth',2, 'DisplayName','Invest');
plot(timesp, cstate, 'LineWidth',2, 'DisplayName','Consumption');
plot(timesp, ystate, 'LineWidth',2, 'DisplayName','Production');
yline(hss, 'DisplayName','Hours SS')
yline(xss, 'DisplayName','Invest SS')
yline(css, 'DisplayName','Consumption SS')
yline(yss, 'DisplayName','Production SS')
yline(kss, 'DisplayName','Captial SS')
legend;
hold off

saveas(figure(3),'simul_5000.jpg')

%% PART B. State-Space Representation and Estimation.
% at = [logkt, logzt, Tlt, Txt, loggt, 1]
% yt = [logyt, logxt, loght, loggt]
% T: 6*6 matrix
% K: 6*4 matrix
% log(kt+1) = A*{log(kt)-log(kss)} + B*St

T = [[A, B, 0];
     [zeros(4,1), P, zeros(4,1)];
     [0, 0, 0, 0, 0, 1]];
 
row_pol = [A, B, 0];

K = [[0, 0, 0, 0];
     Q;
    [0, 0, 0, 0]];

% Log-linearize Resource const, Invest, and Intra
%-----------------------
% Production
syms y k h z
prod = y - (k^theta)*((z*h)^(1-theta));
grad_prod = gradient(prod, [y, k, h, z]);
Px = double(subs(grad_prod, [y, k, h, z], [yss, kss, hss, z_ss]));
Phi_y = [Px(1)*yss; Px(2)*kss; Px(3)*hss; Px(4)*z_ss];
clear y k h z

%------------------------
% Intra-temporal equation
syms h z tl g kp k 
l_int = psi*(((k^theta)*((z*h)^(1-theta))) - grate*kp + (1-delta)*k - g) - ((1-h)*(1-tl)*(1-theta)*((k/h)^theta)*(z^(1-theta)));
grad_int = gradient(l_int, [k, kp, h, z, tl, g]);
Ix = double(subs(grad_int, [k, kp, h, z, tl, g], [kss, kss, hss, z_ss, tl_ss, g_ss]));
Phi_h = [Ix(1)*kss; Ix(2)*kss; Ix(3)*hss; Ix(4)*z_ss; Ix(5)*tl_ss; Ix(6)*g_ss];
clear c h z tl g kp k 

%-- Assortment
% log(yt) = a1*log(zt) + a2*log(kt) + a3*log(ht)
sub_y = [Phi_y(4)/(-Px(1)*yss), Phi_y(2)/(-Px(1)*yss), Phi_y(3)/(-Px(1)*yss)];

% log(ht) = b1*log(kt) + b2*log(kt+1) + b3*log(zt) + b4*log(tl_t) + b5*log(g_ss)
sub_h = [Phi_h(1)/(-Phi_h(3)), Phi_h(2)/(-Phi_h(3)), Phi_h(4)/(-Phi_h(3)), Phi_h(5)/(-Phi_h(3)), Phi_h(6)/(-Phi_h(3))];

% Then,log(yt) = a1*log(zt) + a2*log(kt) + a3*[b1*log(kt) + b2*log(kt+1)...
%                + b3*log(zt) + b4*log(tl_t) + b5*log(g_ss)]
%
%            = (a2+a3*b1)*log(kt) + (a1+a3*b3)*log(zt) + a3*b4*log(tl_t) ...
%                + a3*b5*log(g_ss)+ a3*b2*log(kt+1) 
sub_yh = [(sub_y(2)+sub_y(3)*sub_h(1)), (sub_y(1)+sub_y(3)*sub_h(3)), ...
            sub_y(3)*sub_h(4), sub_y(3)*sub_h(5), sub_y(3)*sub_h(3)];

syms kp k x
invst = x + (1-delta)*k - grate*kp;
grad_inv = gradient(invst, [x, k, kp]);
Ivx = double(subs(grad_inv, [x, k, kp], [xss, kss, kss]));
Phi_x = [Ivx(1)*xss, Ivx(2)*kss, Ivx(3)*kss];

% log(xt) = c1*log(kt) + c2*log(kt+1)
sub_x = [Phi_x(2)/(-Phi_x(1)), Phi_x(3)/(-Phi_x(1))];
        
% 
sub_yh
sub_x
sub_h

C = [[sub_yh(1), sub_yh(2), sub_yh(3), 0, sub_yh(4), 0];
     [ sub_x(1),         0,         0, 0,         0, 0];
     [ sub_h(1),  sub_h(3),  sub_h(4), 0,  sub_h(5), 0];
     [        0,         0,         0, 0,         1, 0]] + [sub_yh(5); sub_x(2); sub_h(2); 0]*row_pol;
 
% Export result of A, B, C
writematrix(T,'A_SSR.txt','Delimiter','tab')
writematrix(K,'B_SSR.txt','Delimiter','tab')
writematrix(C,'C_SSR.txt','Delimiter','tab')
        
%% These are Matrices on State-Space Representation
% X(t+1) = AX(t) + B*eps(t+1)
% Y(t) = C*X(t)
% X(t) = [logkt, logzt, logTlt, logTxt, loggt, 1]
% Y(t) = [logyt, logxt, loght, loggt]

% Re-label it
A = T;
B = K;
S0 = S_bar';
N_iter = 5000;
% Check whether this is correct.
mu = zeros(6, 1, N_iter);
simul = zeros(4, 1, N_iter);
mu(:, 1) = A * [log(1.8)-log(kss); S0'; 1];
simul(:, 1) = C*mu(:, 1);

for i = 2:N_iter
    eps = mvnrnd([0,0,0,0, 0, 0], B*B');
    mu(:, i) = A*mu(:, i-1) + eps';
    simul(:, i) = C*mu(:, i);
end

clear srK srZ srY srH

for i = 1:N_iter
    srK(i) = mu(1, :, i);
    srZ(i) = mu(1, :, i);
    srY(i) = simul(1, :, i);
    srH(i) = simul(3, :, i);
end

timesp = linspace(1, N_iter, N_iter);

figure(4)
plot(timesp, exp(srK+log(kss)), 'LineWidth',2, 'DisplayName','Capital')
title('Simulation of Capital, Shock, Hours, Production with SSR')
xlabel('Time')
ylabel('Level')
hold on
plot(timesp, exp(srZ+log(z_ss)), 'LineWidth',2, 'DisplayName','Shock');
plot(timesp, exp(srY+log(yss)), 'LineWidth',2, 'DisplayName','Production');
plot(timesp, exp(srH+log(hss)), 'LineWidth',2, 'DisplayName','Hours');
yline(hss, 'DisplayName','Hours SS')
yline(z_ss, 'DisplayName','Shock SS')
yline(yss, 'DisplayName','Production SS')
yline(kss, 'DisplayName','Capital SS')
legend;
hold off
saveas(figure(4),'A_test.jpg')

writematrix(mu,'Xt.txt','Delimiter','tab')
writematrix(simul,'Yt.txt','Delimiter','tab')


%% Maximum Likelihood Estimation with Hill-Climbing 

% Initial State Guess
X0 = [mean(mu(1,1,:));
      mean(mu(2,1,:));
      mean(mu(3,1,:));
      mean(mu(4,1,:));
      mean(mu(5,1,:));
      1]; 
S0 = repmat(0.001, 6, 6);    

% Initial parameter guess
startvals = [P; Q];  

% Do Nelson-Mead Algorithm/ fminsearch
kal_fun = @(x) kal_filter(x, simul, X0, S0);
options = optimset('MaxFunEvals', 10000000);
sol = fminsearch(kal_fun, startvals, options);
      
writematrix(sol(1:4,:),'P_MLE.txt','Delimiter','tab')
writematrix(sol(5:8,:),'Q_MLE.txt','Delimiter','tab')

%% Reconstruct the series

% Re-simulation with estimated parameter
temp1 = importdata('A_SSR.txt');
temp2 = importdata('P_MLE.txt');
temp3 = importdata('Q_MLE.txt');

A_hat = [temp1(1,:);
         [zeros(4,1), temp2, zeros(4,1)];
         [zeros(1,5), 1]];
     
B_hat = [zeros(1,4);
         temp3;
         zeros(1,4)];
V_hat = B_hat*B_hat';
C_hat = importdata('C_SSR.txt');

mu_hat = zeros(6, 1, N_iter);
simul_hat = zeros(4, 1, N_iter);

for i = 1:N_iter
    mu_hat(:,:,i) = A_hat*mu(:,:,i) + mvnrnd([0,0,0,0,0,0], V_hat)';
    simul_hat(:,:,i) = C_hat*mu(:,:,i); 
    k_truth(i) = mu(1,1,i);
    k_esti(i) = mu_hat(1,1,i);
    
    x_truth(i) = simul(2,1,i);
    x_esti(i) = simul_hat(2,1,i);
    
    h_truth(i) = simul(3,1,i);
    h_esti(i) = simul_hat(3,1,i);
    
    y_truth(i) = simul(1,1,i);
    y_esti(i) = simul_hat(1,1,i);
end

figure(5)
tiledlayout(2,2);
nexttile
plot(timesp, exp(k_truth+log(kss)), 'LineWidth',1, 'DisplayName','Capital (Truth)')
title('Simulated Capital, Truth versus Estimate')
xlabel('Time')
ylabel('Level')
hold on
plot(timesp, exp(k_esti+log(kss)), 'LineWidth',1, 'DisplayName','Capital ')
yline(kss, 'DisplayName','Capital SS')
legend;
hold off

nexttile
plot(timesp, exp(x_truth+log(xss)), 'LineWidth',1, 'DisplayName','Invest (Truth)')
title('Simulated Invest, Truth versus Estimate')
xlabel('Time')
ylabel('Level')
hold on
plot(timesp, exp(x_esti+log(xss)), 'LineWidth',1, 'DisplayName','Invest')
yline(xss, 'DisplayName','Invest SS')
legend;
hold off

nexttile
plot(timesp, exp(h_truth+log(hss)), 'LineWidth',1, 'DisplayName','Hour (Truth)')
title('Simulated Hours, Truth versus Estimate')
xlabel('Time')
ylabel('Level')
hold on
plot(timesp, exp(h_esti+log(hss)), 'LineWidth',1, 'DisplayName','Hour')
yline(hss, 'DisplayName','Labour SS')
legend;
hold off

nexttile
plot(timesp, exp(y_truth+log(yss)), 'LineWidth',1, 'DisplayName','Production (Truth)')
title('Simulated Production, Truth versus Estimate')
xlabel('Time')
ylabel('Level')
hold on
plot(timesp, exp(y_esti+log(yss)), 'LineWidth',1, 'DisplayName','Production')
yline(yss, 'DisplayName','Production SS')
legend;
hold off
saveas(figure(5),'Truth_vs_estimate.jpg')

toc;
