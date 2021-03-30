function [ass_opt, W_new, state] = EDG_grid(r, w, param)
    beta = param(1); mu = param(2);
    R = r + 1;
    l_low = param(3); l_high = param(4);
    prod_grid = [l_low, l_high];
    
    % Construct asset grid.
    a_max = 6;
    N_ass = 50;
    N_prod = 2;
    a_grid = linspace(0, a_max, N_ass);

    for i = 1:N_ass
        a_grid(i) = (a_grid(i)^2)/a_max;
    end

    C_next = zeros(N_ass, N_prod);
    edg_grid = zeros(N_ass, N_prod);
    C_tilde = zeros(N_ass, N_prod);
    B_pol = zeros(N_ass, N_prod);

    % Construct the First guess for consumption
    for i = 1:N_ass
        for j = 1:N_prod
            C_next(i, j) = R*a_grid(i) + w*prod_grid(j);
        end
    end

    tol = 0.0000001; iter = 0; maxDiff = 100;

    while (maxDiff > tol)
        iter = iter + 1;
        %fprintf('Current iteration is: %.4f \n', iter);
        C_pol = C_next;
        % Construct Endogeneous Grid (Equal World)
        for j = 1:N_prod
            for i = 1:N_ass
                B_pol(i, j) = beta*(R/2)*((C_next(i,1))^-mu + (C_next(i,2))^-mu);
                C_tilde(i, j) = B_pol(i, j)^(-1/mu);
                edg_grid(i, j) = (C_tilde(i, j) + a_grid(i) - w*prod_grid(j))/R;
            end
        end
        
        % Construct Endogeneous Grid (Unequal World)
        % for i = 1:N_ass
        %        B_pol(i, 1) = beta*R*((0.9*(C_next(i,1))^-mu) + ((0.1)*(C_next(i,2))^-mu));
        %        C_tilde(i, 1) = B_pol(i, 1)^(-1/mu);
        %        edg_grid(i, 1) = (C_tilde(i, 1) + a_grid(i) - w*prod_grid(1))/R;
        % end
        
        %for i = 1:N_ass
        %        B_pol(i, 2) = beta*R*(((0.1)*(C_next(i,1))^-mu) + (0.9*(C_next(i,2))^-mu));
        %        C_tilde(i, 2) = B_pol(i, 2)^(-1/mu);
        %        edg_grid(i, 2) = (C_tilde(i, 2) + a_grid(i) - w*prod_grid(2))/R;
        %end    
        
        
        C_next = moveNext(a_grid, edg_grid, C_tilde, prod_grid, N_ass, N_prod, R, w);
        if iter >= 5000
            break
        end
        maxDiff = norm(C_next - C_pol);
    end
    
    ass_opt = zeros(50, 2);

    for i = 1:N_ass
            ass_opt(i, 1) = R*a_grid(i) + w*prod_grid(1) - C_next(i, 1);
            ass_opt(i, 2) = R*a_grid(i) + w*prod_grid(2) - C_next(i, 2);
    end
%    figure(1)
%    plot(a_grid, ass_opt(:, 1), 'LineWidth',2, 'DisplayName','Low');
%    hold on;
%    plot(a_grid, ass_opt(:, 2), 'LineWidth',2, 'DisplayName','High');
%    plot(a_grid, a_grid,'k', 'DisplayName','45 Degree')
%    legend();
%    title('Optimal Asset Policy function')
%    xlabel('level of asset today')
%    ylabel('level of asset tomorrow')
%    hold off

%    figure(2)
%    plot(a_grid, C_next(:, 1), 'LineWidth',2, 'DisplayName','Low');
%    hold on;
%    plot(a_grid, C_next(:, 2), 'LineWidth',2, 'DisplayName','High');
%    legend();
%    title('Optimal Consumption Policy function')
%    xlabel('level of asset today')
%    ylabel('level of consumption today')
%    hold off

    %% Construct Stationary distribution of asset
    % Monte Carlo for the invariant distribution
    A = [ass_opt(:, 2)', ass_opt(:, 1)'];
    % Last index
    M = 50;
    
    %Initial distribution for Capital holding
    W = 6*rand(10000, 1);
    W_new = zeros(size(W));
    [N, ] = histcounts(W, linspace(0,6,100));
    N = N'/10000;
    N_new = zeros(99,1);

    % generate initial state (L or H)
    state = rand(10000, 1);
    state = (state>0.5); % define state = 1 for High
    maxPdiff = 100;

    while (maxPdiff > 0.01)
        %fprintf('Current iteration is: %.8f \n', maxPdiff);
        
        %N_new = N;
        H = rand(10000,1);
        H = (H<=0.5);
        L = rand(10000,1);
        L = (L<=0.5);
        
         % Indicate state from the initial distribution
         for i = 1:10000
            if state(i)==1
                if H(i)==1
                    state(i)=1;
                else
                    state(i)=0;
                end
            else 
                if L(i)==1
                    state(i)=0;
                else
                    state(i)=1;
                end
            end
        end

        for i=1:10000
            %compute next period wealth if the state is H today 
            if state(i)==1
                if W(i)< a_grid(1)
                    W_new(i)=0;
                else
                    for j=1:M-1
                        if W(i) >= a_grid(j) && W(i) < a_grid(j+1)
                           W_new(i)=(a_grid(j+1)-W(i))/(a_grid(j+1)-a_grid(j))*A(j)...
                               +(W(i) - a_grid(j))/(a_grid(j+1)-a_grid(j))*A(j+1);
                            continue
                        end 
                    end
                    if W(i)>=a_grid(M)
                        W_new(i)=(W(i)-a_grid(M-1))/(a_grid(M)-a_grid(M-1))*A(M);
                    end
                end
            else
                
            %compute next period wealth if the state is L today
                if W(i) < a_grid(1)
                    W_new(i) = 0;
                else
                    for j=1:M-1
                        if W(i)>=a_grid(j) && W(i)<a_grid(j+1)
                            W_new(i)=(a_grid(j+1)-W(i))/(a_grid(j+1)-a_grid(j))*A(M+j)...
                                +(W(i) - a_grid(j))/(a_grid(j+1)-a_grid(j))*A(M+j+1);
                            continue
                        end
                    end
                    if W(i) >= a_grid(M)
                        W_new(i) = (W(i)-a_grid(M-1))/(a_grid(M)-a_grid(M-1))*A(2*M);
                    end
                end
            end
        end
        % update the histogram count
        [N_new, ] = histcounts(W_new,linspace(0,6,100));
        N_new = N_new'/10000;
        maxPdiff = norm(N_new - N);
        N = N_new;
    end
end

