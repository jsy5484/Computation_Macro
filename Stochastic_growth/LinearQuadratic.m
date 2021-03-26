clear all
tic;

% Parameterization
theta = 0.35; delta = 0.0464; gamma_z = 0.016; gamma_n = 0.015;
beta = 0.9722; beta_hat = beta*(1+gamma_n); sigma = 0.5; rho = 0.2; psi = 2.24;
grate = (1+gamma_z)*(1+gamma_n);
param = [theta, delta, gamma_z, gamma_n, beta_hat, psi];

% Find a Steady State
SS = SS(param);
kss = SS(1);
hss = SS(2);
lss = 1 - hss;
zss = 0;
css = ((kss^param(1))*((exp(0)*hss)^(1-param(1)))) - (1+param(4))*(1+param(3))*kss + (1-param(2))*kss;

% Apply Kydland-Prescott with subroutine to find h*. 
syms k kp z 
f = @(k, kp, z) log(((k^theta)*((exp(z)*h_find(k, kp, z, psi, theta, grate, delta))^(1-theta)))-grate*kp+(1-delta)*k) ...
     +psi*log(1-h_find(k, kp, z, psi, theta, grate, delta));

r_bar = f(kss, kss, zss);
r_bar = double(r_bar);
z_bar = [kss, zss, 1, kss];
z_bar = z_bar';

eps = 0.001;
y = [f(kss, kss, zss), f(kss+eps, kss, zss)];
x = [kss, kss+eps];
grad1 = diff(y)/diff(x);
clear y x

y = [f(kss, kss, zss), f(kss, kss, zss+eps)];
x = [zss, zss+eps];
grad2 = diff(y)/diff(x);
clear y x

y = [f(kss, kss, zss), f(kss, kss+eps, zss)];
x = [kss, kss+eps];
grad3 = diff(y)/diff(x);
clear y x

% Construct Gradient
grad = [grad1 grad2 0 grad3];
grad = grad';

% Construct Hessian Matrix
% 1st Row
y = [f(kss, kss, zss), f(kss+eps, kss, zss), f(kss-eps, kss, zss)];
hess11 = (y(2) - 2*y(1)+y(3))/eps^2;
clear y
y = [f(kss, kss, zss), f(kss+eps, kss, zss), f(kss, kss, zss+eps), f(kss+eps, kss, zss+eps)];
hess12 = (-((y(2)-y(1))/eps) + ((y(4)-y(3))/eps))/eps;
clear y
y = [f(kss, kss, zss), f(kss+eps, kss, zss), f(kss, kss+eps, zss), f(kss+eps, kss+eps, zss)];
hess14 = (-((y(2)-y(1))/eps) + ((y(4)-y(3))/eps))/eps;
clear y

% 2nd Row
y = [f(kss, kss, zss), f(kss, kss, zss+eps), f(kss+eps, kss, zss), f(kss+eps, kss, zss+eps)];
hess21 = (-((y(2)-y(1))/eps) +((y(4)-y(3))/eps))/eps;
clear y
y = [f(kss, kss, zss), f(kss, kss, zss+eps), f(kss, kss, zss-eps)];
hess22 = (y(2) - 2*y(1)+y(3))/eps^2;
clear y
y = [f(kss, kss, zss), f(kss, kss, zss+eps), f(kss, kss+eps, zss), f(kss, kss+eps, zss+eps)];
hess24 = (-((y(2)-y(1))/eps) + ((y(4)-y(3))/eps))/eps;
clear y

% 4nd Row
y = [f(kss, kss, zss), f(kss, kss+eps, zss), f(kss+eps, kss, zss), f(kss+eps, kss+eps, zss)];
hess41 = (-((y(2)-y(1))/eps) + ((y(4)-y(3))/eps))/eps;
clear y
y = [f(kss, kss, zss), f(kss, kss+eps, zss), f(kss, kss, zss+eps), f(kss, kss+eps, zss+eps)];
hess42 = (-((y(2)-y(1))/eps) + ((y(4)-y(3))/eps))/eps;
clear y
y = [f(kss, kss, zss), f(kss, kss+eps, zss), f(kss, kss-eps, zss)];
hess44 = (y(2) - 2*y(1)+y(3))/eps^2;
clear y

hess = [hess11 hess12 0 hess14
        hess21 hess22 0 hess24
        0      0      0 0     
        hess41 hess42 0 hess44];
    
clear a b c d h11 h12 h13 h14 h21 h22 h23 h24 h31 h32 h33 h34 h41 h42 h43 h44
ev = [0,0,1,0];
ev = ev';

in1 = r_bar - grad.'*z_bar + (1/2)*z_bar.'*hess*z_bar;
A = z_bar.'*hess;
B = hess*z_bar;
in2 = grad*ev.' - ev*A - B*ev.' + ev*grad.';

M = ev*in1*ev.' + (1/2)*in2 + (1/2)*hess;

Q = [[M(1,1), M(1,2), M(1,3)];
     [M(2,1), M(2,2), M(2,3)];
     [M(3,1), M(3,2), M(3,3)]];

R = M(4,4);

W = [M(1,4);
     M(2,4);
     M(3,4)];
 
A = [[0,0,0];
     [0,rho,0];
     [0,0,1]];
 
B = [1;
     0;
     0];
 
C = [0; 
     1;
     0];

% Transform into undiscounted problem
A_tld = sqrt(beta_hat)*(A-B*inv(R)*W.');
B_tld = sqrt(beta_hat)*B;
Q_tld = Q - W*inv(R)*W.';

cH11 = inv(A_tld);
cH12 = inv(A_tld)*B_tld*inv(R)*B_tld.';
cH21 = Q_tld * inv(A_tld);
cH22 = Q_tld*inv(A_tld)*B_tld*inv(R)*B_tld.' + A_tld.';

H = [[cH11, cH12];
     [cH21, cH22]];

[eigVec, eigVal] = eig(H); 
v = nonzeros(eigVal');
eigVal = reshape(v, 6,1);

Gam = [[eigVal(4), 0, 0, 0, 0, 0];
       [0, eigVal(5), 0, 0, 0, 0];
       [0, 0, eigVal(6), 0, 0, 0];
       [0, 0, 0, eigVal(1), 0, 0];
       [0, 0, 0, 0, eigVal(2), 0];
       [0, 0, 0, 0, 0, eigVal(3)]];
   
V = [eigVec(:,4),eigVec(:,5),eigVec(:,6),eigVec(:,1),eigVec(:,2),eigVec(:,3)];
V_inv = inv(V);

W22 = [[V_inv(4,4), V_inv(4,5), V_inv(4,6)];
       [V_inv(5,4), V_inv(5,5), V_inv(5,6)];
       [V_inv(6,4), V_inv(6,5), V_inv(6,6)]];
   
W21 = [[V_inv(4,1), V_inv(4,2), V_inv(4,3)];
       [V_inv(5,1), V_inv(5,2), V_inv(5,3)];
       [V_inv(6,1), V_inv(6,2), V_inv(6,3)]];

P = -inv(W22)*W21;
inner = R + B_tld.'*P*B_tld;
F = inv(inner)*B_tld.'*P*A_tld + inv(R)*W.';


% Constructing Policy function
k_grid = linspace(0.5*kss, 1.5*kss, 100);
k_pol_L = zeros(100, 1);
k_pol_H = zeros(100, 1);
for i = 1:100
    k_pol_L(i) = -F*[k_grid(i);
                             -0.5; 
                              1];
end

for i = 1:100
    k_pol_H(i) = -F*[k_grid(i);
                       0.5; 
                       1];
end


figure(1)
plot(k_grid, k_pol_L, 'LineWidth',2, 'DisplayName','Low');
hold on;
plot(k_grid, k_pol_H, 'LineWidth',2, 'DisplayName','High');
legend();
xlim([0.5*kss 1.5*kss])
title('Optimal Capital Policy Function from Linear-Quadratic_(c)')
xlabel('level of Capital Today')
ylabel('level of Capital Tomorrow')
hold off
saveas(figure(1),'LQ_C_capital.jpg');

toc;

