function Xt = gen_IRF(A, B, T)
% This function generate an impulse response using A and B.

Xt  = zeros(size(A,1), T+1);
eps = zeros(size(B,2), T+1);
eps(1,2) = 0.1;

for t = 2:T+1
	Xt(:,t) = A*Xt(:,t-1) + B*eps(:,t);
end

Xt = Xt(:,2:51);
end
