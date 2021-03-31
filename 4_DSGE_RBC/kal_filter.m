function loglhood = kal_filter(x, simul, X0, S0)

N_iter = length(simul);
temp1 = importdata('A_SSR.txt');


A = [temp1(1,:);
     [0,      x(1, 1:4), 0];
     [0,      x(2, 1:4), 0];
     [0,      x(3, 1:4), 0];
     [0,      x(4, 1:4), 0];
     [zeros(1,5), 1]];
 
B = [ zeros(1,4);
      x(5:8, 1:4);
       zeros(1,4)];
Bt = B*B';

C = importdata('C_SSR.txt');

X = zeros(6, 1, N_iter);
X_cond = zeros(6, 1, N_iter);
S = zeros(6, 6, N_iter);
S_cond = zeros(6, 6, N_iter);
m_err = zeros(4, 1, N_iter); 
F_mat = zeros(4, 4, N_iter);
logll = zeros(N_iter, 1);

X_esti = X0;
S_esti = S0;

for i = 1:N_iter
    % Time update
    X_next = A*X_esti;
    S_next = A*S_esti*A' + Bt;
   
    X_cond(:,1,i) = X_next;
    S_cond(:,:,i) = S_next;
    
    % Measurement 
    m_err(:,1,i) = ([simul(1,i); simul(2,i); simul(3,i); simul(4,i)] - C*X_next);
    F_mat(:,:,i) = C*S_next*C';
    X_esti = X_next + S_next*C'*inv(C*S_next*C')*m_err(:,1,i);
    S_esti = S_next - S_next*C'*inv(C*S_next*C')*C*S_next;
   
    X(:,1,i) = X_esti;
    S(:,:,i) = S_esti;
    
    logll(i) = (1/2)*(log(det(F_mat(:,:,i)))+m_err(:,1,i)'*inv(F_mat(:,:,i))*m_err(:,1,i));
end


loglhood = (N_iter/2)*log(2*pi) + sum(logll);
fprintf('Current log-likelihood is: %.4f \n', loglhood);

end

