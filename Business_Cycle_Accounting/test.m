
% Parameterization
theta = 0.3333334; delta = 0.05; gamma_z = 0.016; gamma_n = 0.015;
beta = 0.9722; beta_hat = beta*(1+gamma_n); psi = 2.5;

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
yss = (kss^theta)*((z_ss*hss)^(1-theta));


%% STEP 2. % Log-linearize the Intratemporal equation
syms h z tl g kp k 
c = (k^theta)*((z*h)^(1-theta)) - grate*kp + (1-delta)*k - g;
w = (1-theta)*((k/h)^theta)*(z^(1-theta));

l_int = psi*(((k^theta)*((z*h)^(1-theta))) - grate*kp + (1-delta)*k - g) - ((1-h)*(1-tl)*(1-theta)*((k/h)^theta)*(z^(1-theta)));
%l_int = (1-tl)*w*(1-h)-psi*c;
grad_int = gradient(l_int, [k, kp, h, z, tl, g]);
Ix = double(subs(grad_int, [k, kp, h, z, tl, g], [kss, kss, hss, z_ss, tl_ss, g_ss]));

coeff = [Ix(1)*kss;
         Ix(2)*kss;
         Ix(3)*hss;
         Ix(4)*z_ss;
         Ix(5)*tl_ss;
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

l_Eul = -beta_hat*c*(rp + (1-delta)*(1+txp)) + cp*(1+tx)*grate;


grad_Eul = gradient(l_Eul, [k, kp, kpp, h, hp, z, tx, g, zp, txp, gp]);
Ex = double(subs(grad_Eul, [k, kp, kpp, h, hp, z, tx, g, zp, txp, gp], [kss, kss, kss, hss, hss, z_ss, tx_ss, g_ss, z_ss, tx_ss, g_ss]));


% 11 arguments
coeff_Eul = [Ex(1)*kss;
             Ex(2)*kss;
             Ex(3)*kss;
             Ex(4)*hss;
             Ex(5)*hss;
             Ex(6)*z_ss;
             Ex(7)*tx_ss;
             Ex(8)*g_ss;
             Ex(9)*z_ss;
             Ex(10)*tx_ss;
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
%% Simulation
% Simulate for 1000 periods
horiz = 800;
timesp = linspace(1, horiz, horiz);

St = zeros(4,horiz);
% Setup starting point
St(:,1) = zeros(4,1);

for i = 2:horiz
    eps = mvnrnd([0 0 0 0], eye(4));
    St(:,i) =  P*St(:,i-1)  + Q*eps';
end

zhat  = exp(St(1,:))*z_ss; tlhat = exp(St(2,:))*tl_ss;
txhat = exp(St(3,:))*tx_ss; ghat  = exp(St(4,:))*g_ss;

plot(1:1:800, zhat); hold on; plot(1:1:800, tlhat);
plot(1:1:800, txhat); plot(1:1:800, ghat); hold off; 

llkhat(1) = 0;

for i = 1:horiz
    llkhat(i+1) = A*llkhat(i)  +  B*St(:,i);
    llhhat(i)   = C2*llkhat(i) + D2*St(:,i);
end

khat = exp(llkhat)*kss;
hhat = exp(llhhat)*hss;

for i = 1:horiz
    xhat(i) = grate*khat(i+1) - (1-delta)*khat(i);
end

for i = 1:horiz
    yhat(i) = (khat(i)^theta)*((zhat(i)*hhat(i))^(1-theta));
    chat(i) = yhat(i) - xhat(i) - ghat(i);
end
llkhat(end)=[];


plot(1:1:800, (chat - css)./css, 'LineWidth',1.5, 'DisplayName','Consumption'); hold on;
plot(1:1:800, (xhat - xss)./xss, 'LineWidth',1.5,'DisplayName','Investment'); 
plot(1:1:800, (yhat - yss)./yss, 'LineWidth',1.5, 'DisplayName','Output');
plot(1:1:800, llhhat, 'LineWidth',1.5, 'DisplayName','Hours'); legend;