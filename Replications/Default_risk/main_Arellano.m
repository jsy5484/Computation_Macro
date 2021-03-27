% This is a Matlab program for solving Arellano's Income fluctuation problem

clear all; clc; tic;
%% Parameterization
beta = 0.953; sigma = 2; r = 0.017; rho = 0.945; eta = 0.025; theta = 0.282;
ny = 21; nB = 400; tol = 1e-8; 

%% Construct grids
% Grid for Bonds
bmax = 0.5; bmin = -0.5; db = (bmax - bmin)/(nB - 1);
b=bmin:db:bmax; b(find(b==min(b(b>=0)))) = 0; b=b(:);
bzero = find(b==0);   

% Discretization
[Z,Zprob1] = tauchenHussey(ny, 0, rho, eta, eta); ygrid = exp(Z);
ydefgrid = min(0.969 * mean(ygrid), ygrid);
Zprob = kron(1, Zprob1);

% Initial Guess
q0 = 1/(1 + r)*ones(nB*ny, nB);
Vc = ones(nB*ny, 1); Vd = ones(nB*ny, 1);  
Vdc = reshape(Vc, nB, ny*1); Vdc = Vdc(bzero,:);      
Vdc = ones(nB, 1)*Vdc;         
Vdc = reshape(Vdc, nB*ny*1, 1); 

y = kron(ones(1,1), ygrid).*kron(1, ones(ny, 1));
y = kron(y, ones(nB, 1));

ydef = kron(ones(1,1), ydefgrid).*kron(1, ones(ny,1));
ydef = kron(ydef, ones(nB,1)); 

%% Main Iteration for a price functional
diff_r = 1; tol_r = 1e-6; iter = 0;
while diff_r>tol_r
    iter = iter + 1;
    % Saving Function
    S = (kron(1, ones(nB*ny, 1))*b').*q0-(kron(ones(ny*1, 1),b)*ones(1, nB));
    
    % Current consumption when not in default 
    c = y*ones(1,nB) - S;
    % Current consumption when in default 
    cdefault = ydef*ones(1, nB);

    % Calculate utility first
    if (sigma ~= 1)                            
        u = (((c).^(1 - sigma)))./(1 - sigma);
        u(find(c<=0)) = NaN;                   
        udefault = (((cdefault).^(1 - sigma)))./(1 - sigma);
        udefault(find(cdefault<=eps)) = NaN;
    elseif (sigma == 1)                         
        u = log(max(c, realmin));
        u(find(c<=realmin)) = NaN;             
        udefault = log(max(cdefault, eps));
        udefault(find(cdefault<=eps)) = NaN;
    end
    
    % Inner iteration for VFI
    diff = 1000; tolV = min(max(diff_r, tol_r), 1e-6); 
    while diff>tolV
        % Evaluate expected value based on each state
        EVc = Zprob*(reshape(Vc,nB,ny)');   % Expected value when in the Contract
        EVc = kron(EVc,ones(nB,1)); 

        EVd = Zprob*(reshape(Vd,nB,ny*1)'); % Expected value when in Default
        EVd = kron(EVd,ones(nB,1));

        EVdc = Zprob*(reshape(Vdc,nB,ny*1)'); % Expected value when enter in Contract after had defaulted
        EVdc = kron(EVdc,ones(nB,1));

        [Vd1, policybad] = max(udefault + beta.*theta.*EVdc + beta.*(1-theta).*EVd,[],2);
        [Vc1, policygood] = max(u + beta.*EVc,[],2);
        default = Vd1>Vc1 | isnan(Vc1) == 1;
        Vc1(find(Vd1>Vc1 | isnan(Vc1) == 1)) = Vd1(find(Vd1>Vc1| isnan(Vc1) == 1));

        % Find new Values
        Vdc = reshape(Vc1,nB,ny*1);
        Vdc = Vdc(bzero,:);
        Vdc = ones(nB,1)*Vdc;
        Vdc = reshape(Vdc,nB*ny*1,1);


        % Update
        diff = max([max(max(abs(Vc - Vc1))), max(max(abs(Vd - Vd1)))]); 
        Vc = Vc1;
        Vd = Vd1;
    end 
    
    fprintf('Current Iteration is: %.d \n', iter);
    fprintf('Convergence of Value is obtained. Now updating a Price functional \n');
    
    % Find a new price and update
    % Expected default
    Edef = Zprob*(reshape(default, nB, ny*1)'); 
    ind = reshape(default, nB, ny*1)';
    ind = sum(ind)==ny*1; 

    for j = 1:nB
        if ind(j)==1
            Edef(:,j)=1;
        end
    end

    Edef = kron(Edef,ones(nB, 1)); 
    q1 = (1/(1+r))*((1 - Edef));
    q1 = max(q1, 0);

    diff_r = max(max(abs(q1-q0))); 
    q0 = q1;
end

V = max(Vc, Vd);
Sprime = (kron(1,ones(nB*ny,1))*b').*q0-(kron(ones(ny*1,1),b)*ones(1,nB));

%% Find converged policy functions
policy=(1-default).*policygood + default.*bzero;

for i = 1:nB*ny
    ct(i,1) = c(i,policy(i));
    qt(i,1) = q0(i,policy(i));
end

toc;
%% Graph
high = mean(ygrid).*1.05; low = mean(ygrid).*0.95;
% Find the place where high and low are contained in the ygrid. Actually is
% the closest value in the ygrid, since it cannot contain the same value
i_low = find(ygrid==max(ygrid(find(abs(low - ygrid)==min(abs(low - ygrid))))));
i_high = find(ygrid==max(ygrid(find(abs(high - ygrid)==min(abs(high - ygrid))))));
yhigh = i_high*nB; ylow = i_low*nB;

% Probability of Default 
figure1 = figure(1);
titulo=sprintf('Probability of Default');
mesh(b,ygrid,abs(reshape(default,nB,ny))');
axis([min(b),0,min(ygrid),max(ygrid),0,1]);
view(0,90);
xlabel('B')
ylabel('y')
colormap(flipud(bone));
title(titulo);

% Bond price schedules 
figure(2);
plot(b, q0(yhigh,:),'LineWidth',2, 'DisplayName', 'y high');
hold on; grid on;
plot(b, q0(ylow,:),'LineWidth',2, 'DisplayName', 'y low');
xlim([-0.35 0]); title('Bond prices schedules - q(B,y)')
xlabel('B'); ylabel('q(B,y)'); legend;
hold off;

% Equilibrium interest rate 
sovr = 1./q0 -1; sovr(sovr==Inf) = nan; 
sovrhigh = sovr(yhigh,:); sovrlow = sovr(ylow,:);

figure(3)
plot(b, sovrhigh,'LineWidth',2, 'DisplayName', 'y high');
hold on; grid on;
plot(b, sovrlow,'LineWidth',2, 'DisplayName', 'y low');
xlim([-0.13 0.02]); ylim([0 0.20]);
title('Equilibrium interest rate - 1/q(B,y)')
xlabel('B'); legend; ylabel('1/q(B,y)');
hold off;

% Savings Function 
default1 = 0;
policy1 = (1 - default1).*policygood + default1.*bzero;
% Bprime is represented by
bprime = (kron(1,ones(nB*ny,1))*b').*q0;

% Calculate bprime_t(state)
for i=1:nB*ny
    bprimet(i,1)=bprime(i,policy1(i));
end

figure(4);
plot(b, bprimet(yhigh:((nB-1)+yhigh)),'LineWidth',2, 'DisplayName', 'y high');
hold on; grid on;
plot(b, bprimet(ylow :((nB-1)+ylow)),'LineWidth',2, 'DisplayName', 'y high');
xlim([-0.2 0.15]); ylim([-0.15 0.10]); xlabel('B'); ylabel('Savings');
title('Saving Function'); legend;
hold off;

% Value Function 
V_low = V(ylow:((nB - 1) + ylow));
V_high = V(yhigh:((nB - 1) + yhigh));

figure(5);
plot(b, V_high,'LineWidth',2, 'DisplayName', 'y high');
hold on; grid on;
plot(b, V_low,'LineWidth',2, 'DisplayName', 'y low');
xlim([-0.4 0.15]); xlabel('B'); ylabel('v^0(B,y)')
legend; title('Value Function - v^0(B,y)');
hold off


%Simulation();


 
