function retrn = simult(A, B, Xss, horiz)
% This function simulates of Xt using A and B, and return K, w, r.

global Params; rng(100);
eps = normrnd(0, Params.sigma^2, horiz,1);

X_diff      = zeros(size(A,1),1,horiz+1);
Xt          = zeros(size(A,1),1,horiz+1);
Xt(:,:,1)   = Xss;

for i = 1:horiz
    X_diff(:,:,i+1) = A*X_diff(:,:,i) + B*eps(i);
    Xt(:,:,i+1)     = Xss + X_diff(:,:,i+1);
end

for i = 1:horiz
    K_agg(i) = Params.knotDistrK'*Xt(61:1059,i) + Params.knotDistrK'*Xt(1060:2058,i);
    R(i) = intr((K_agg(i)/0.9300), exp(Xt(2059,1,i)));
    W(i) = cwage((K_agg(i)/0.9300), exp(Xt(2059,1,i)));
end

retrn = [K_agg; R; W];
end

