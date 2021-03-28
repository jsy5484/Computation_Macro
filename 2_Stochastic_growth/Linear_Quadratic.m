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
css = ((kss^param(1))*((exp(0)*hss)^(1-param(1)))) - (1+param(4))*(1+param(3))*kss + (1-param(2))*kss;

syms k kp z h
f = log(((k^param(1))*((exp(z)*h)^(1-param(1)))) - (1+param(4))*(1+param(3))*kp + (1-param(2))*k)+param(6)*log(1-h);

r_bar = subs(f, [k, kp, z, h], [kss, kss, 0, hss]);
r_bar = double(r_bar);
z_bar = [kss, 0, 1, kss, hss];
z_bar = z_bar';

a = diff(f, k); b = diff(f, z); c = diff(f, kp); d = diff(f, h);
grad = [double(subs(a, [k, kp, z, h], [kss, kss, 0, hss])), double(subs(b, [k, kp, z, h], [kss, kss, 0, hss])), 0, double(subs(c, [k, kp, z, h], [kss, kss, 0, hss])), double(subs(d, [k, kp, z, h], [kss, kss, 0, hss]))];
grad = grad';

% Construct Hessian Matrix
h11 = diff(a, k); h12 = diff(a, z); h13 = diff(a, kp); h14 = diff(a, h);
h21 = diff(b, k); h22 = diff(b, z); h23 = diff(b, kp); h24 = diff(b, h);
h31 = diff(c, k); h32 = diff(c, z); h33 = diff(c, kp); h34 = diff(c, h);
h41 = diff(d, k); h42 = diff(d, z); h43 = diff(d, kp); h44 = diff(d, h);

hess = [[double(subs(h11, [k, kp, z, h], [kss, kss, 0, hss])), double(subs(h12, [k, kp, z, h], [kss, kss, 0, hss])), 0, double(subs(h13, [k, kp, z, h], [kss, kss, 0, hss])), double(subs(h14, [k, kp, z, h], [kss, kss, 0, hss]))];
        [double(subs(h21, [k, kp, z, h], [kss, kss, 0, hss])), double(subs(h22, [k, kp, z, h], [kss, kss, 0, hss])), 0, double(subs(h23, [k, kp, z, h], [kss, kss, 0, hss])), double(subs(h24, [k, kp, z, h], [kss, kss, 0, hss]))];
        [0, 0, 0, 0, 0];
        [double(subs(h31, [k, kp, z, h], [kss, kss, 0, hss])), double(subs(h32, [k, kp, z, h], [kss, kss, 0, hss])), 0, double(subs(h33, [k, kp, z, h], [kss, kss, 0, hss])), double(subs(h34, [k, kp, z, h], [kss, kss, 0, hss]))];
        [double(subs(h41, [k, kp, z, h], [kss, kss, 0, hss])), double(subs(h42, [k, kp, z, h], [kss, kss, 0, hss])), 0, double(subs(h43, [k, kp, z, h], [kss, kss, 0, hss])), double(subs(h44, [k, kp, z, h], [kss, kss, 0, hss]))]];

clear a b c d h11 h12 h13 h14 h21 h22 h23 h24 h31 h32 h33 h34 h41 h42 h43 h44
ev = [0,0,1,0,0];
ev = ev';

in1 = r_bar - grad.'*z_bar + (1/2)*z_bar.'*hess*z_bar;
A = z_bar.'*hess;
B = hess*z_bar;
in2 = grad*ev.' - ev*A - B*ev.' + ev*grad.';

M = ev*in1*ev.' + (1/2)*in2 + (1/2)*hess;

Q = [[M(1,1), M(1,2), M(1,3)];
     [M(2,1), M(2,2), M(2,3)];
     [M(3,1), M(3,2), M(3,3)]];

R = [[M(4,4),M(4,5)];
     [M(5,4),M(5,5)]];

W = [[M(1,4),M(1,5)];
     [M(2,4),M(2,5)];
     [M(3,4),M(3,5)]];
 
A = [[0,0,0];
     [0,rho,0];
     [0,0,1]];
 
B = [[1,0];
     [0,0];
     [0,0]];
 
C = [0; 
     1;
     0];

A_tld = sqrt(beta_hat)*(A-B*inv(R)*W.');

B_tld = sqrt(beta_hat)*B;

Q_tld = Q - W*inv(R)*W.';


P0 = [[0,0,0];
      [0,0,0];
      [0,0,0]];
  
tol = 1.0e-10;

[Pn, Fn] = Ricatti(P0, A_tld, B_tld, R, Q_tld, tol);

F = Fn + inv(R)*W.';
P = Pn;

% Constructing Policy function
k_grid = linspace(0.5*kss, 1.5*kss, 100);
pol_L = zeros(2, 1, 100);
for i = 1:100
    pol_L(:,:,i) = -F*[k_grid(i);
                             -0.5; 
                              1];
end

k_pol_L =zeros(1, 100);
h_pol_L =zeros(1, 100);

for i = 1:100
    k_pol_L(1, i) = pol_L(1,:,i);
    h_pol_L(1, i) = pol_L(2,:,i);

end

X_H = linspace(0, 7, 100);
pol_H = zeros(2, 1, 100);
for i = 1:100
    pol_H(:,:,i) = -F*[k_grid(i);
                       0.5; 
                       1];
end
k_pol_H =zeros(1, 100);
h_pol_H =zeros(1, 100);
for i = 1:100
    k_pol_H(1, i) = pol_H(1,:,i);
    h_pol_H(1, i) = pol_H(2,:,i);
end

figure(1)
plot(k_grid, k_pol_L, 'LineWidth',2, 'DisplayName','Low');
hold on;
plot(k_grid, k_pol_H, 'LineWidth',2, 'DisplayName','High');
legend();
xlim([0.5*kss 1.5*kss])
title('Optimal Capital Policy Function from Linear-Quadratic_(b)')
xlabel('level of Capital Today')
ylabel('level of Capital Tomorrow')
hold off
saveas(figure(1),'LQ_B_capital.jpg');

figure(2)
plot(k_grid, h_pol_L, 'LineWidth',2, 'DisplayName','Low');
hold on;
plot(k_grid, h_pol_H, 'LineWidth',2, 'DisplayName','High');
legend();
xlim([0.5*kss 1.5*kss])
title('Optimal Labour Policy Function from Linear-Quadratic_(b)')
xlabel('level of Capital Today')
ylabel('level of Labour')
hold off
saveas(figure(2),'LQ_B_hrs.jpg');

toc;
