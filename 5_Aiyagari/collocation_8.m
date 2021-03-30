clear all
clc; tic; 
% parameterization
beta = 0.96; mu = 3; w = 1.0; R = 1.02; chi = 5000; 

% Grid over saving
a_max = 6;
ngrid_a = 8;
a_grid = linspace(0.01, a_max, ngrid_a);

% Grid over productivity
ngrid_z = 2;
l_low = 0.2; l_high = 0.9;
z_grid = [l_low, l_high];
param = [mu R w chi l_low l_high];

% Construct Basis Function
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

% Initial Guess for theta with vectorization
% # of nodes = 16 (2*8)
theta = zeros(1, 2*ngrid_a);
for i = 1:ngrid_a
    theta(1, i) = 0.9*(R*a_grid(i)+w*l_low);
    theta(1, i+ngrid_a) = 0.8*(R*a_grid(i)+w*l_high);
end

%% Main Loop
iter = 0; max_iter = 50; diff = 100; tol = 1e-2;

while diff >tol && iter < max_iter
    iter = iter + 1;
    fprintf('Current iteration is : %d  \n', iter);
    clear temp Residt Resid agrid_val agrid_valt a_tt a_tp
    a_pol_col = zeros(1, ngrid_a*2);
    a_pol_colt = zeros(1, ngrid_a*2);
    Resid = zeros(1, 2*ngrid_a);
  
    % Construct the policy function with initial guess
    a_pol_L= theta(1)*psi_1_L + theta(2)*psi_2_L + theta(3)*psi_3_L + theta(4)*psi_4_L +theta(5)*psi_5_L + theta(6)*psi_6_L+ ...
             theta(7)*psi_7_L + theta(8)*psi_8_L;
    
    a_pol_H= theta(9)*psi_1_H + theta(10)*psi_2_H + theta(11)*psi_3_H + theta(12)*psi_4_H +theta(13)*psi_5_H + theta(14)*psi_6_H+ ...
             theta(15)*psi_7_H + theta(16)*psi_8_H;
         
    % Evaluate the policy and Residual function on collocation points
    for i = 1:ngrid_a
        % Evaluate the policy function on collocation
        a_pol_col(i)         = double(subs(a_pol_L, a, a_grid(i)));
        a_pol_col(i+ngrid_a) = double(subs(a_pol_H, a, a_grid(i)));
        
        % Evaluate the Residual function on collocation 
        Resid(i)         = ((R*a_grid(i) + w*l_low - a_pol_col(i))^-mu)         -(beta*R/2)*expected(a_pol_col(i),         a_pol_L, a_pol_H, param);
        Resid(i+ngrid_a) = ((R*a_grid(i) + w*l_high- a_pol_col(i+ngrid_a))^-mu) -(beta*R/2)*expected(a_pol_col(i+ngrid_a), a_pol_L, a_pol_H, param);
    end
    
    % We need value of a_t such that making all of the value of Resid to 0.
    
    % Now updating with Netwon
    % Calculate Numerical Jacobian + updating rule 
    eps = 0.01;
    temp = zeros(1, 2*ngrid_a);
    Residt = zeros(2*ngrid_a, 2*ngrid_a);
    thetat = zeros(1, 2*ngrid_a);

    for i = 1:2*ngrid_a
        temp = zeros(1, 2*ngrid_a);
        temp(i) = eps;
        thetat = theta + temp;  
        %temp
        a_pol_Lt= thetat(1)*psi_1_L + thetat(2)*psi_2_L + thetat(3)*psi_3_L + thetat(4)*psi_4_L +thetat(5)*psi_5_L + thetat(6)*psi_6_L+ ...
                  thetat(7)*psi_7_L + thetat(8)*psi_8_L;
              
        a_pol_Ht= thetat(9)*psi_1_H + thetat(10)*psi_2_H + thetat(11)*psi_3_H + thetat(12)*psi_4_H  + thetat(13)*psi_5_H  + thetat(14)*psi_6_H+ ...
                  thetat(15)*psi_7_H + thetat(16)*psi_8_H;
              

        for j = 1:ngrid_a
            a_pol_colt(j) = double(subs(a_pol_Lt, a, a_grid(j)));
            a_pol_colt(j+ngrid_a) = double(subs(a_pol_Ht, a, a_grid(j)));
            
            Residt(j,i)         =((R*a_grid(j) + w*l_low  - a_pol_colt(j))^-mu        ) - (beta*R/2)*expected(a_pol_colt(j),         a_pol_Lt, a_pol_Ht, param);
            Residt(j+ngrid_a,i) =((R*a_grid(j) + w*l_high - a_pol_colt(j+ngrid_a))^-mu) - (beta*R/2)*expected(a_pol_colt(j+ngrid_a), a_pol_Lt, a_pol_Ht, param);
            %Residt
        end
    end
    
    J = zeros(2*ngrid_a, 2*ngrid_a);

    Resid = Resid';
    J = [(Residt(:,1 )-Resid)/eps, (Residt(:,2 )-Resid)/eps, (Residt(:,3 )-Resid)/eps, (Residt(:,4 )-Resid)/eps, (Residt(:,5 )-Resid)/eps,...
         (Residt(:,6 )-Resid)/eps, (Residt(:,7 )-Resid)/eps, (Residt(:,8 )-Resid)/eps, (Residt(:,9 )-Resid)/eps, (Residt(:,10)-Resid)/eps,...
         (Residt(:,11)-Resid)/eps, (Residt(:,12)-Resid)/eps, (Residt(:,13)-Resid)/eps, (Residt(:,14)-Resid)/eps, (Residt(:,15)-Resid)/eps,...
         (Residt(:,16)-Resid)/eps];
    a_tp = theta' - inv(J)*Resid;
    diff = norm(a_tp-theta');
    fprintf('Current Difference is : %.6f \n', diff);

    theta = a_tp;
    theta = theta';

end

writematrix(theta,'FEM.txt','Delimiter','tab')

toc;



                                                  
       