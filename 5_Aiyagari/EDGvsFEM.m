clear all
clc

% Parameterization
beta = 0.96; R = 1.02; w = 1.0; mu = 3.0; l_low = 0.2; l_high = 0.9;
param = [beta mu];
alpha = 0.33; r = 0.02;

% Finite Element Method with Collocation
a_max = 6;
ngrid_acol = 8;
a_grid = linspace(0.01, a_max, ngrid_acol);
FEM = importdata('FEM.txt');

% Defining Piecewise Basis function
% Low
syms a 
psi_1_L = piecewise(a_grid(1)<=a<=a_grid(2), ((a_grid(2)  - a)/(a_grid(2)-a_grid(1))), 0);
psi_2_L = piecewise(a_grid(1)<=a<=a_grid(2), ((a -  a_grid(1))/(a_grid(2)-a_grid(1))), a_grid(2)<=a<=a_grid(3), ( (a_grid(3) - a)/(a_grid(3)-a_grid(2)) ), 0); 
psi_3_L = piecewise(a_grid(2)<=a<=a_grid(3), ((a -  a_grid(2))/(a_grid(3)-a_grid(2))), a_grid(3)<=a<=a_grid(4), ( (a_grid(4) - a)/(a_grid(4)-a_grid(3)) ), 0);
psi_4_L = piecewise(a_grid(3)<=a<=a_grid(4), ((a -  a_grid(3))/(a_grid(4)-a_grid(3))), a_grid(4)<=a<=a_grid(5), ( (a_grid(5) - a)/(a_grid(5)-a_grid(4)) ), 0);             
psi_5_L = piecewise(a_grid(4)<=a<=a_grid(5), ((a -  a_grid(4))/(a_grid(5)-a_grid(4))), a_grid(5)<=a<=a_grid(6), ( (a_grid(6) - a)/(a_grid(6)-a_grid(5)) ), 0);
psi_6_L = piecewise(a_grid(5)<=a<=a_grid(6), ((a -  a_grid(5))/(a_grid(6)-a_grid(5))), a_grid(6)<=a<=a_grid(7), ( (a_grid(7) - a)/(a_grid(7)-a_grid(6)) ), 0);
psi_7_L = piecewise(a_grid(6)<=a<=a_grid(7), ((a -  a_grid(6))/(a_grid(7)-a_grid(6))), a_grid(7)<=a<=a_grid(8), ( (a_grid(8) - a)/(a_grid(8)-a_grid(7)) ), 0);
psi_8_L = piecewise(a_grid(7)<=a<=a_grid(8), ((a -  a_grid(7))/(a_grid(8)-a_grid(7))), 0);
% High
psi_1_H = piecewise(a_grid(1)<=a<=a_grid(2), ((a_grid(2)-a)/(a_grid(2)-a_grid(1))), 0);
psi_2_H = piecewise(a_grid(1)<=a<=a_grid(2), ((a -  a_grid(1))/(a_grid(2)-a_grid(1))), a_grid(2)<=a<=a_grid(3), ( (a_grid(3) - a)/(a_grid(3)-a_grid(2)) ), 0); 
psi_3_H = piecewise(a_grid(2)<=a<=a_grid(3), ((a -  a_grid(2))/(a_grid(3)-a_grid(2))), a_grid(3)<=a<=a_grid(4), ( (a_grid(4) - a)/(a_grid(4)-a_grid(3)) ), 0);
psi_4_H = piecewise(a_grid(3)<=a<=a_grid(4), ((a -  a_grid(3))/(a_grid(4)-a_grid(3))), a_grid(4)<=a<=a_grid(5), ( (a_grid(5) - a)/(a_grid(5)-a_grid(4)) ), 0);             
psi_5_H = piecewise(a_grid(4)<=a<=a_grid(5), ((a -  a_grid(4))/(a_grid(5)-a_grid(4))), a_grid(5)<=a<=a_grid(6), ( (a_grid(6) - a)/(a_grid(6)-a_grid(5)) ), 0);
psi_6_H = piecewise(a_grid(5)<=a<=a_grid(6), ((a -  a_grid(5))/(a_grid(6)-a_grid(5))), a_grid(6)<=a<=a_grid(7), ( (a_grid(7) - a)/(a_grid(7)-a_grid(6)) ), 0);
psi_7_H = piecewise(a_grid(6)<=a<=a_grid(7), ((a -  a_grid(6))/(a_grid(7)-a_grid(6))), a_grid(7)<=a<=a_grid(8), ( (a_grid(8) - a)/(a_grid(8)-a_grid(7)) ), 0);
psi_8_H = piecewise(a_grid(7)<=a<=a_grid(8), ((a -  a_grid(7))/(a_grid(8)-a_grid(7))), 0);

a_pol_L_star = FEM(1)*psi_1_L + FEM(2)*psi_2_L + FEM(3)*psi_3_L + FEM(4)*psi_4_L + FEM(5)*psi_5_L + FEM(6)*psi_6_L +...
               FEM(7)*psi_7_L + FEM(8)*psi_8_L;
           
a_pol_H_star = FEM(9)*psi_1_H + FEM(10)*psi_2_H + FEM(11)*psi_3_H + FEM(12)*psi_4_H + FEM(13)*psi_5_H  + FEM(14)*psi_6_H +...
               FEM(15)*psi_7_H + FEM(16)*psi_8_H;


figure(1)
fplot(a_pol_L_star,'LineWidth', 2, 'DisplayName','Low');
title('Optimal Saving function with FEM')
ylim([-0.0185 7])
xlim([0 6.3])
legend('Location','northwest');

xlabel('level of asset today')
ylabel('level of asset tomorrow')
grid on
hold on
fplot(a_pol_H_star, 'LineWidth', 2, 'DisplayName','High');
plot(a_grid, FEM(1:8), '^', 'MarkerSize',7,'MarkerFaceColor','black', 'DisplayName','coefficient of basis (Low)');
plot(a_grid, FEM(9:16), 'square', 'MarkerSize',7,'MarkerFaceColor','black', 'DisplayName','coefficient of basis (High)');
hold off
saveas(figure(1), 'opt_sav_FEM.jpg')


% Adjust grid for EGM
a_max = 6;
ngrid_a = 50;
a_grid_edm = linspace(0.01, a_max, ngrid_a);

for i = 1:ngrid_a
    a_grid_edm(i) = (a_grid_edm(i)^2)/a_max;
end

% Solve the problem with Endogenous Grid Method
[ass_opt, W, state] = EDG_grid(r,w, param);

figure(2)
plot(a_grid_edm, ass_opt(:, 1), 'LineWidth', 2, 'DisplayName','Low (EGM)');
grid on;
hold on;
plot(a_grid_edm, ass_opt(:, 2), 'LineWidth', 2, 'DisplayName','High (EGM)');
plot(a_grid_edm, a_grid_edm,'k', 'DisplayName','45 Degree')
ylim([-0.019, 7])
title('Comparison of Optimal saving function, FEM vs EGM')
xlabel('level of asset today')
ylabel('level of asset tomorrow')
plot(a_grid, FEM(1:8), '^', 'MarkerSize',7,'MarkerFaceColor','black', 'DisplayName','Low (FEM)');
plot(a_grid, FEM(9:16), 'square', 'MarkerSize',7,'MarkerFaceColor','black', 'DisplayName','High (FEM)');
legend('Location','northwest');
saveas(figure(2), 'FEMvsEGM.jpg')
hold off

