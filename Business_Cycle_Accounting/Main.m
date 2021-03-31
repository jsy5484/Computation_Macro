clear all; clc
tic;

global Params; initparam;

% Model Matched Sbar : Satisfying Long-run Mean
z_ss = exp(0.002454975848088);
tl_ss = 0.116886124345528;
tx_ss = 0.133788352975186;
g_ss = 0.1528;


% Function find_Sb is a rountine that find z_bar, tl_bar, tx_bar such that
% kss, xss, yss, hss matches with long-run mean.
% Number above is obatined by running this function. 

%X0 = [0; 0.1; 0.9];
%findS = @(X) find_Sb(X, ihat, hhat, yhat, ghat);
%options = optimset('MaxFunEvals', inf);
%soln = fminsearch(findS, X0, options);

% Initial Guess
P       = .995*eye(4);
Q = zeros(4); V = Q*Q';
Sbar    = [log(z_ss); tl_ss; tx_ss; log(g_ss)];
P0      = (eye(4)-P)*Sbar; 
param   = [Params.gn;Params.gz;Params.beta;Params.delta;Params.psi;...
           Params.sigma;Params.theta; Params.qadj];

%% STEP 1: Given Sbar, find Steady-State.
% Find Steady-State of Economy       
SS = BCA_findSS(Sbar);
ks = SS(1); cs = SS(2); hs = SS(3); xs = SS(4); ys = SS(5); zs = SS(6);
tauls = SS(7); tauxs = SS(8); gs = SS(9);

%% STEP 2: Solve this economy : Find A, B and C matrix.
X0 = [log(ks);log(zs);tauls;tauxs;log(gs);1];
Y0 = [log(ys);log(xs);log(hs);log(gs)];
[A, B, C] = BCA_solve(SS, param, P, Q, P0, X0, Y0); 

% Export result of A, B, C
writematrix(A,'A_SSR.txt','Delimiter','tab')
writematrix(B,'B_SSR.txt','Delimiter','tab')
writematrix(C,'C_SSR.txt','Delimiter','tab')

%% STEP 3: Import Data: Base Year 1990.
Yhat = BCA_data(Params.gz);
yhat = Yhat(:,1); ihat = Yhat(:,2); hhat = Yhat(:,3); ghat = Yhat(:,4);
ytilde = log(yhat) - log(mean(yhat));
itilde = log(ihat) - log(mean(ihat));
htilde = log(hhat) - log(mean(hhat));
gtilde = log(ghat) - log(mean(ghat));

%Ytilde = [ytilde, itilde, htilde, gtilde]';
Ytilde = [log(yhat), log(ihat), log(hhat), log(ghat)]';

%% STEP 4. Kalman Filtering and implement MLE.
% Initial parameter guess
startvals = [P; 0.1*eye(4)];   
S0 = repmat(0.001, 6, 6);    
 
% Do Nelson-Mead Algorithm/ fminsearch
kal_fun = @(x) kal_filter(x, Ytilde, X0, S0);
options = optimset('MaxFunEvals', inf);
sol = fminsearch(kal_fun, startvals, options);

P = sol(1:4, 1:4); Q = sol(5:8, 1:4);
startvals = [P; Q];
P0      = (eye(4)-P)*Sbar; 

%writematrix(P,'P.txt','Delimiter','tab')
%writematrix(Q,'Q.txt','Delimiter','tab')
P = importdata('P.txt')
Q = importdata('Q.txt')
% Solve again using the estimated P and Q.
X0 = [log(ks); log(zs); tauls; tauxs; log(gs); 1];
Y0 = [log(ys); log(xs); log(hs); log(gs)];
[A, B, C] = BCA_solve(SS, param, P, Q, P0, X0, Y0); 
%% STEP 5. Extract Wedges
[extr_St, lkhat] = BCA_extract(C, Sbar, Yhat);
lzhat  = extr_St(:,1);
ltlhat = extr_St(:,2);
ltxhat = extr_St(:,3);
lghat  = extr_St(:,4);
Tlhat  = extr_St(:,5);
Zhat   = extr_St(:,6);

horiz = length(extr_St(:,1));
lyhat = Ytilde(1,:)'; lxhat = Ytilde(2,:)'; lhhat = Ytilde(3,:)';
%lghat = Ytilde(4,:)';

%% Last Step : Generate Each Wedge Alone Economy.
% Year starting from 2008. (idx 117 = 2008 Q1)
Xt0 = [lkhat(117:end), lzhat(117:end), ltlhat(117:end), ltxhat(117:end), lghat(117:end), ones(40,1)];
YM0 = Xt0*C';

% Re-estimate the decision rule, C-matrix by fixing expectation to control
% the wedges turning on and off.
Y0 = 1;
s0 = [lzhat(Y0); ltlhat(Y0); ltxhat(Y0); lghat(Y0)];
[~, ~, C0] = sol_fixexp(Sbar, P, Q, s0, [0;0;0;0]); 
[~, ~, C1] = sol_fixexp(Sbar, P, Q, s0, [1;0;0;0]); 
[~, ~, C2] = sol_fixexp(Sbar, P, Q, s0, [0;1;0;0]);
[~, ~, C3] = sol_fixexp(Sbar, P, Q, s0, [0;0;1;0]);
[~, ~, C4] = sol_fixexp(Sbar, P, Q, s0, [0;0;0;1]);
temp = ones(40,1);

% No Weges / Only one : Efficiency / Labour / Investment / Government Exp
YMn      = (Xt0 - temp*Xt0(Y0,:))*C0' + temp*YM0(Y0,:);
YMz      = (Xt0 - temp*Xt0(Y0,:))*(C1 - C0)' + temp*YM0(Y0,:);
YMl      = (Xt0 - temp*Xt0(Y0,:))*(C2 - C0)' + temp*YM0(Y0,:);
YMx      = (Xt0 - temp*Xt0(Y0,:))*(C3 - C0)' + temp*YM0(Y0,:);
YMg      = (Xt0 - temp*Xt0(Y0,:))*(C4 - C0)' + temp*YM0(Y0,:);

% All Weges / All but : Efficiency / Labour / Investment / Government Exp
YMall    = (Xt0 - temp*Xt0(Y0,:))*C' + temp*YM0(Y0,:);
YMnoz    = (Xt0 - temp*Xt0(Y0,:))*(C2 + C3 + C4 -2*C0)' + temp*YM0(Y0,:);
YMnol    = (Xt0 - temp*Xt0(Y0,:))*(C1 + C3 + C4 -2*C0)' + temp*YM0(Y0,:);
YMnox    = (Xt0 - temp*Xt0(Y0,:))*(C1 + C2 + C4 -2*C0)' + temp*YM0(Y0,:);
Y_actual = yhat(117:end);

toc;
%% Draw Graphes
figure(1)
plot(1:1:156, 100*exp(lyhat-lyhat(1)), 'LineWidth',1.5, 'DisplayName','Actual Output');
hold on; title('U.S. output and three measured wedges','FontSize',13 ); 
ylabel('Index : 1979 Q1 = 100')
plot(1:1:156, 100*(Zhat/Zhat(1)).^(1-Params.theta), 'LineWidth',1.5, 'DisplayName','Extracted Efficienty Wedge');
plot(1:1:156, 100*(1-Tlhat)/(1-Tlhat(1)), 'LineWidth',1.5, 'DisplayName','Extracted Labour Wedge');
plot(1:1:156, 100*((1+ltxhat(Y0))*ones(horiz,1)./(1+ltxhat)),'LineWidth',1.5, 'DisplayName','Extracted Investment Wedge');
legend('Location','Northwest'); grid on;
xticks([1 17 33 49 65 81 97 113 129 145 156])
xticklabels({'1979-Q1', '1983-Q1', '1987-Q1', '1991-Q1', '1995-Q1', '1999-Q1',...
      '2003-Q1', '2007-Q1', '20011-Q1', '2015-Q1', '2017-Q3'}); xtickangle(45)
  
figure(2)
plot(1:1:40, 100*exp(Y_actual-Y_actual(1)), 'LineWidth',1.5, 'DisplayName','Actual GDP');
hold on; title('Data and prediction with One Wedges in Great Depression','FontSize',13 ); 
ylabel('Index : 2008 Q1 = 100')
plot(1:1:40, 100*exp(YMz(:,1)-YMz(Y0,1)), 'LineWidth',1.5, 'DisplayName','Efficiency Wedge Only');
plot(1:1:40, 100*exp(YMl(:,1)-YMl(Y0,1)),  'LineWidth',1.5, 'DisplayName','Labour Wedge Only');
plot(1:1:40, 100*exp(YMx(:,1)-YMx(Y0,1)), 'LineWidth',1.5, 'DisplayName','Investment Wedge Only');
legend('Location','Best'); xlim([1 40]); grid on; xticks([1 3 11 19 27 35 40])
xticklabels({'2008-Q1', '2009-Q1', '2011-Q1', '2013-Q1', '2015-Q1', '2017-Q1', '2017-Q3'})
xtickangle(45)
hold off


figure(3)
plot(1:1:40, 100*exp(Y_actual-Y_actual(1)), 'LineWidth',1.5, 'DisplayName','Actual GDP');
hold on; title('Data and prediction with All but One wedge in Great Depression','FontSize',13); 
ylabel('Index : 2008 Q1 = 100')
plot(1:1:40, 100*exp(YMnoz(:,1)-YMnoz(Y0,1)), 'LineWidth',1.5, 'DisplayName','All but No Efficiency Wedge');
plot(1:1:40, 100*exp(YMnol(:,1)-YMnol(Y0,1)),  'LineWidth',1.5, 'DisplayName','All but No Labour Wedge');
plot(1:1:40, 100*exp(YMnox(:,1)-YMnox(Y0,1)), 'LineWidth',1.5, 'DisplayName','All but No Investment Wedge');
legend('Location','Best'); xlim([1 40]); grid on; xticks([1 3 11 19 27 35 40])
xticklabels({'2008-Q1', '2009-Q1', '2011-Q1', '2013-Q1', '2015-Q1', '2017-Q1', '2017-Q3'})
xtickangle(45)
hold off
  
