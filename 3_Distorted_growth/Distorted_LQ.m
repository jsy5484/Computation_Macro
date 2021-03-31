clear all
clc
tic;

% Solving tax-distorted economy with LQ confusing. 

% Parameterization
theta = 0.35; delta = 0.0464; gamma_z = 0.016; gamma_n = 0.015;
beta = 0.9722; beta_hat = beta*(1+gamma_n); psi = 2.24;
param = [theta, delta, gamma_z, gamma_n, beta_hat, psi];
grate = (1+gamma_z)*(1+gamma_n);

P =[[ 0.9887,  0.0307, -0.0089,  -0.0407];
    [-0.0012,  1.0011, -0.0275,   0.0175];
    [-0.0045,  0.0449,  0.9675,  -0.0426];
    [ 0.0063,  0.0017,  0.0016,   0.9945]];

P0 = [0.0140; 0.0008; 0.0129; -0.0137];

S_bar = inv(eye(4)-P)*P0;
z_ss = exp(S_bar(1));
tl_ss = exp(S_bar(2));
tx_ss = exp(S_bar(3));
g_ss = exp(S_bar(4));

%% STEP 1.
% Find a Steady State
% Kss, Css, Hss obatained analytically.
kl_rto = ( ((1+tx_ss)*(grate-beta_hat*(1-delta))) / (beta_hat*theta*(z_ss^(1-theta))))^(1/(theta-1));
zt1 = (((kl_rto)^(theta-1))*(z_ss^(1-theta))-grate+1-delta);
zt2 = (1/psi)*((1-tl_ss)*(1-theta)*((kl_rto)^theta)*(z_ss^(1-theta)));
zt3 = zt2/kl_rto;

kss = (zt2+g_ss)/(zt1+zt3);
css = zt1*kss - g_ss;
hss = (1/kl_rto)*kss;
xss = grate*kss - (1-delta)*kss;

% Adjust for logged wedges
z_ss = S_bar(1);
tl_ss = S_bar(2);
tx_ss = S_bar(3);
g_ss = S_bar(4);

syms k z tl tx g AK AH AX kp h 
% Express Return function
r = theta*(AK^(theta-1))*((exp(z)*AH)^(1-theta));
w = (1-theta)*(AK^theta)*((exp(z)*AH)^(-theta))*exp(z);

f = log(r*k + (1-exp(tl))*w*h - (1+exp(tx))*((grate*kp - (1-delta)*k)) + exp(tl)*w*AH + exp(tx)*AX - exp(g)) + psi*log(1-h);
r_bar = subs(f, [k, z, tl, tx, g, AK, AH, AX, kp, h], [kss, z_ss, tl_ss, tx_ss, g_ss, kss, hss, xss, kss, hss]);
r_bar = double(r_bar);
z_bar = [kss, 1, z_ss, tl_ss, tx_ss, g_ss, kss, xss, hss, kss, hss];
% Order   k,   z,     tl,    tx,    g,   1,   K,   X,  H, || k,   h
% This is 12 elements vector
z_bar = z_bar';

a = diff(f, k);
% We have 1
b = diff(f, z);
c = diff(f, tl);
d = diff(f, tx);
e = diff(f, g);
%---------------------
f_ = diff(f, AK);
g_ = diff(f, AX);
h_ = diff(f, AH);
%---------------------
i_ = diff(f, kp);
j_ = diff(f, h);


var1 = [k, kp, h, z, tl, tx, g, AK, AX, AH];
ss_pt = [kss, kss, hss, z_ss, tl_ss, tx_ss, g_ss, kss, xss, hss];

grad = [double(subs(a, var1, ss_pt)),...
        0,...
        double(subs(b, var1, ss_pt)),...
        double(subs(c, var1, ss_pt)),...
        double(subs(d, var1, ss_pt)),...
        double(subs(e, var1, ss_pt)),...
        double(subs(f_, var1, ss_pt)),...
        double(subs(g_, var1, ss_pt)),...
        double(subs(h_, var1, ss_pt)),...
        double(subs(i_, var1, ss_pt)),...
        double(subs(j_, var1, ss_pt))];
grad = grad';

% This is a very inefficient way to get Hessian, but I wanted to check the matrix visually, because of the location of 1 matters.
% For the rest of code, even though it's a bit time-consuming, I made everything visually.
% Construct Hessian Matrix
h_a1 = diff(a, k); h_a2 = diff(a, z); h_a3 = diff(a, tl); h_a4 = diff(a, tx);
h_a5 = diff(a, g); h_a6 = diff(a, AK); h_a7 = diff(a, AX); h_a8 = diff(a, AH);
h_a9 = diff(a, kp); h_a10 = diff(a, h);

h_b1 = diff(b, k); h_b2 = diff(b, z); h_b3 = diff(b, tl); h_b4 = diff(b, tx);
h_b5 = diff(b, g); h_b6 = diff(b, AK); h_b7 = diff(b, AX); h_b8 = diff(b, AH);
h_b9 = diff(b, kp); h_b10 = diff(b, h);

h_c1 = diff(c, k); h_c2 = diff(c, z); h_c3 = diff(c, tl); h_c4 = diff(c, tx);
h_c5 = diff(c, g); h_c6 = diff(c, AK); h_c7 = diff(c, AX); h_c8 = diff(c, AH);
h_c9 = diff(c, kp); h_c10 = diff(c, h);

h_d1 = diff(d, k); h_d2 = diff(d, z); h_d3 = diff(d, tl); h_d4 = diff(d, tx);
h_d5 = diff(d, g); h_d6 = diff(d, AK); h_d7 = diff(d, AX); h_d8 = diff(d, AH);
h_d9 = diff(d, kp); h_d10 = diff(d, h);

h_e1 = diff(e, k); h_e2 = diff(e, z); h_e3 = diff(e, tl); h_e4 = diff(e, tx);
h_e5 = diff(e, g); h_e6 = diff(e, AK); h_e7 = diff(e, AX); h_e8 = diff(e, AH);
h_e9 = diff(e, kp); h_e10 = diff(e, h);

h_f1 = diff(f_, k); h_f2 = diff(f_, z); h_f3 = diff(f_, tl); h_f4 = diff(f_, tx);
h_f5 = diff(f_, g); h_f6 = diff(f_, AK); h_f7 = diff(f_, AX); h_f8 = diff(f_, AH);
h_f9 = diff(f_, kp); h_f10 = diff(f_, h);

h_g1 = diff(g_, k); h_g2 = diff(g_, z); h_g3 = diff(g_, tl); h_g4 = diff(g_, tx);
h_g5 = diff(g_, g); h_g6 = diff(g_, AK); h_g7 = diff(g_, AX); h_g8 = diff(g_, AH);
h_g9 = diff(g_, kp); h_g10 = diff(g_, h);

h_h1 = diff(h_, k); h_h2 = diff(h_, z); h_h3 = diff(h_, tl); h_h4 = diff(h_, tx);
h_h5 = diff(h_, g); h_h6 = diff(h_, AK); h_h7 = diff(h_, AX); h_h8 = diff(h_, AH);
h_h9 = diff(h_, kp); h_h10 = diff(h_, h);

h_i1 = diff(i_, k); h_i2 = diff(i_, z); h_i3 = diff(i_, tl); h_i4 = diff(i_, tx);
h_i5 = diff(i_, g); h_i6 = diff(i_, AK); h_i7 = diff(i_, AX); h_i8 = diff(i_, AH);
h_i9 = diff(i_, kp); h_i10 = diff(i_, h);

h_j1 = diff(j_, k); h_j2 = diff(j_, z); h_j3 = diff(j_, tl); h_j4 = diff(j_, tx);
h_j5 = diff(j_, g); h_j6 = diff(j_, AK); h_j7 = diff(j_, AX); h_j8 = diff(j_, AH);
h_j9 = diff(j_, kp); h_j10 = diff(j_, h);

hess = [[double(subs(h_a1, var1, ss_pt)),   0, double(subs(h_a2, var1, ss_pt)), double(subs(h_a3, var1, ss_pt)), double(subs(h_a4, var1, ss_pt)), double(subs(h_a5, var1, ss_pt)), double(subs(h_a6, var1, ss_pt)), double(subs(h_a7, var1, ss_pt)), double(subs(h_a8, var1, ss_pt)), double(subs(h_a9, var1, ss_pt)), double(subs(h_a10, var1, ss_pt))];
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];
        [double(subs(h_b1, var1, ss_pt)),   0, double(subs(h_b2, var1, ss_pt)), double(subs(h_b3, var1, ss_pt)), double(subs(h_b4, var1, ss_pt)), double(subs(h_b5, var1, ss_pt)), double(subs(h_b6, var1, ss_pt)), double(subs(h_b7, var1, ss_pt)), double(subs(h_b8, var1, ss_pt)), double(subs(h_b9, var1, ss_pt)), double(subs(h_b10, var1, ss_pt))];
        [double(subs(h_c1, var1, ss_pt)),   0, double(subs(h_c2, var1, ss_pt)), double(subs(h_c3, var1, ss_pt)), double(subs(h_c4, var1, ss_pt)), double(subs(h_c5, var1, ss_pt)), double(subs(h_c6, var1, ss_pt)), double(subs(h_c7, var1, ss_pt)), double(subs(h_c8, var1, ss_pt)), double(subs(h_c9, var1, ss_pt)), double(subs(h_c10, var1, ss_pt))];
        [double(subs(h_d1, var1, ss_pt)),   0, double(subs(h_d2, var1, ss_pt)), double(subs(h_d3, var1, ss_pt)), double(subs(h_d4, var1, ss_pt)), double(subs(h_d5, var1, ss_pt)), double(subs(h_d6, var1, ss_pt)), double(subs(h_d7, var1, ss_pt)), double(subs(h_d8, var1, ss_pt)), double(subs(h_d9, var1, ss_pt)), double(subs(h_d10, var1, ss_pt))];
        [double(subs(h_e1, var1, ss_pt)),   0, double(subs(h_e2, var1, ss_pt)), double(subs(h_e3, var1, ss_pt)), double(subs(h_e4, var1, ss_pt)), double(subs(h_e5, var1, ss_pt)), double(subs(h_e6, var1, ss_pt)), double(subs(h_e7, var1, ss_pt)), double(subs(h_e8, var1, ss_pt)), double(subs(h_e9, var1, ss_pt)), double(subs(h_e10, var1, ss_pt))];
        [double(subs(h_f1, var1, ss_pt)),   0, double(subs(h_f2, var1, ss_pt)), double(subs(h_f3, var1, ss_pt)), double(subs(h_f4, var1, ss_pt)), double(subs(h_f5, var1, ss_pt)), double(subs(h_f6, var1, ss_pt)), double(subs(h_f7, var1, ss_pt)), double(subs(h_f8, var1, ss_pt)), double(subs(h_f9, var1, ss_pt)), double(subs(h_f10, var1, ss_pt))];
        [double(subs(h_g1, var1, ss_pt)),   0, double(subs(h_g2, var1, ss_pt)), double(subs(h_g3, var1, ss_pt)), double(subs(h_g4, var1, ss_pt)), double(subs(h_g5, var1, ss_pt)), double(subs(h_g6, var1, ss_pt)), double(subs(h_g7, var1, ss_pt)), double(subs(h_g8, var1, ss_pt)), double(subs(h_g9, var1, ss_pt)), double(subs(h_g10, var1, ss_pt))];
        [double(subs(h_h1, var1, ss_pt)),   0, double(subs(h_h2, var1, ss_pt)), double(subs(h_h3, var1, ss_pt)), double(subs(h_h4, var1, ss_pt)), double(subs(h_h5, var1, ss_pt)), double(subs(h_h6, var1, ss_pt)), double(subs(h_h7, var1, ss_pt)), double(subs(h_h8, var1, ss_pt)), double(subs(h_h9, var1, ss_pt)), double(subs(h_h10, var1, ss_pt))];
        [double(subs(h_i1, var1, ss_pt)),   0, double(subs(h_i2, var1, ss_pt)), double(subs(h_i3, var1, ss_pt)), double(subs(h_i4, var1, ss_pt)), double(subs(h_i5, var1, ss_pt)), double(subs(h_i6, var1, ss_pt)), double(subs(h_i7, var1, ss_pt)), double(subs(h_i8, var1, ss_pt)), double(subs(h_i9, var1, ss_pt)), double(subs(h_i10, var1, ss_pt))];
        [double(subs(h_j1, var1, ss_pt)),   0, double(subs(h_j2, var1, ss_pt)), double(subs(h_j3, var1, ss_pt)), double(subs(h_j4, var1, ss_pt)), double(subs(h_j5, var1, ss_pt)), double(subs(h_j6, var1, ss_pt)), double(subs(h_j7, var1, ss_pt)), double(subs(h_j8, var1, ss_pt)), double(subs(h_j9, var1, ss_pt)), double(subs(h_j10, var1, ss_pt))]];

clear a b c d h11 h12 h13 h14 h21 h22 h23 h24 h31 h32 h33 h34 h41 h42 h43 h44
ev = [0, 1, 0, 0, 0, 0, 0, 0, 0, 0 ,0];
ev = ev';

in1 = r_bar - grad.'*z_bar + (1/2)*z_bar.'*hess*z_bar;
A = z_bar.'*hess;
B = hess*z_bar;
in2 = grad*ev.' - ev*A - B*ev.' + ev*grad.';

M = ev*in1*ev.' + (1/2)*in2 + (1/2)*hess;

% Now we get LQ ingredient
Q = [M(1,1:9);
     M(2,1:9);
     M(3,1:9);
     M(4,1:9);
     M(5,1:9);
     M(6,1:9);
     M(7,1:9);
     M(8,1:9);
     M(9,1:9);];
  
R = [M(10,10:11);
     M(11,10:11)];

W = [M(1,10:11);
     M(2,10:11);
     M(3,10:11);
     M(4,10:11);
     M(5,10:11);
     M(6,10:11);
     M(7,10:11);
     M(8,10:11);
     M(9,10:11);];
 
A = [[0,        0,       0,      0,       0,       0,   0,0,0];
     [0,        1,       0,      0,       0,       0,   0,0,0];
     [0,   0.0140,  0.9887, 0.0307, -0.0089, -0.0407,   0,0,0];
     [0,   0.0008, -0.0012, 1.0011, -0.0275,  0.0175,   0,0,0];
     [0,   0.0129, -0.0045, 0.0449,  0.9675, -0.0426,   0,0,0];
     [0,  -0.0137,  0.0063, 0.0017,  0.0016,  0.9945,   0,0,0];
     [0,        0,       0,      0,       0,        0,  0,0,0];
     [0,        0,       0,      0,       0,        0,  0,0,0];
     [0,        0,       0,      0,       0,        0,  0,0,0]];
  
B = [[1,0];
     [0,0];
     [0,0];
     [0,0];
     [0,0];
     [0,0];
     [0,0];
     [0,0];
     [0,0];];
 
C = [[0,0,       0,       0,       0,       0,  0,0];
     [0,0,       0,       0,       0,       0,  0,0];
     [0,0,  0.0077,  0.0024, -0.0041,  0.0003,  0,0];
     [0,0,  0.0024,  0.0043,  0.0023,  0.0153,  0,0];
     [0,0, -0.0041,  0.0023,  0.0088,  0.0121,  0,0];
     [0,0,  0.0003,  0.0153,  0.0121,  0.0139,  0,0];
     [0,0,0,0,0,0,0,0];
     [0,0,0,0,0,0,0,0];
     [0,0,0,0,0,0,0,0];
     [0,0,0,0,0,0,0,0]];

 
 
%% STEP 2.
% Map this into undiscounted problem.
A_tld = sqrt(beta_hat)*(A-B*inv(R)*W.');
B_tld = sqrt(beta_hat)*B;
Q_tld = Q - W*inv(R)*W.';

%% STEP 3. Contruct ingredient we will use
Q_tld_y = [Q_tld(1,1:6);
           Q_tld(2,1:6);
           Q_tld(3,1:6);
           Q_tld(4,1:6);
           Q_tld(5,1:6);
           Q_tld(6,1:6);];       
       
Q_tld_z = [Q_tld(1,7:9);
           Q_tld(2,7:9);
           Q_tld(3,7:9);
           Q_tld(4,7:9);
           Q_tld(5,7:9);
           Q_tld(6,7:9);];
 
A_tld_y = [A_tld(1,1:6);
           A_tld(2,1:6);
           A_tld(3,1:6);
           A_tld(4,1:6);
           A_tld(5,1:6);
           A_tld(6,1:6);]; 
  
A_z_tld = [A_tld(1,7:9);
           A_tld(2,7:9);
           A_tld(3,7:9);
           A_tld(4,7:9);
           A_tld(5,7:9);
           A_tld(6,7:9);];
       
A_y = [A(1,1:6);
       A(2,1:6);
       A(3,1:6);
       A(4,1:6);
       A(5,1:6);
       A(6,1:6);];


B_tld_y = [B_tld(1,1:2);
           B_tld(2,1:2);
           B_tld(3,1:2);
           B_tld(4,1:2);
           B_tld(5,1:2);
           B_tld(6,1:2);];

B_y = [B(1,1:2);
       B(2,1:2);
       B(3,1:2);
       B(4,1:2);
       B(5,1:2);
       B(6,1:2);];

W_y = [[W(1,1),W(1,2)];
       [W(2,1),W(2,2)];
       [W(3,1),W(3,2)];
       [W(4,1),W(4,2)];
       [W(5,1),W(5,2)];
       [W(6,1),W(6,2)];];
   
W_z = [[W(7,1),W(7,2)];
       [W(8,1),W(8,2)];
       [W(9,1),W(9,2)];];

%% Step 4. Using the Market Clearing condition now.
OMG = [[1, 0, 0, 0, 0, 0];
       [-(1-delta), 0, 0, 0, 0, 0];
       [0, 0, 0, 0, 0, 0];];
   
PSI = [[0,0];
       [grate,0];
       [0,1];];
   
OMG_tld = inv((eye(3) + PSI*inv(R)*W_z.'))*(OMG - PSI*inv(R)*W_y.');
PSI_tld = inv((eye(3) + PSI*inv(R)*W_z.'))*PSI;

A_hat = A_tld_y + A_z_tld*OMG_tld;
Q_hat = Q_tld_y + Q_tld_z*OMG_tld;
B_hat = B_tld_y + A_z_tld*PSI_tld;

A_bar = A_tld_y - B_tld_y*inv(R)*PSI_tld.'*Q_tld_z.';

%% Constructing H matrix
H = [[inv(A_hat), inv(A_hat)*B_hat*inv(R)*B_tld_y.'];
     [Q_hat*inv(A_hat), Q_hat*inv(A_hat)*B_hat*inv(R)*B_tld_y.'+A_bar.']];
  
[eigVec, eigVal] = eig(H); 
v = nonzeros(eigVal');
eigVal = reshape(v, 12,1);

Gam = [[eigVal(7), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];
       [0, eigVal(8), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];
       [0, 0, eigVal(9), 0, 0, 0, 0, 0, 0, 0, 0, 0];
       [0, 0, 0, eigVal(10), 0, 0, 0, 0, 0, 0, 0, 0];
       [0, 0, 0, 0, eigVal(11), 0, 0, 0, 0, 0, 0, 0];
       [0, 0, 0, 0, 0, eigVal(12), 0, 0, 0, 0, 0, 0];
       [0, 0, 0, 0, 0, 0, eigVal(1), 0, 0, 0, 0, 0];
       [0, 0, 0, 0, 0, 0, 0, eigVal(2), 0, 0, 0, 0];
       [0, 0, 0, 0, 0, 0, 0, 0, eigVal(3), 0, 0, 0];
       [0, 0, 0, 0, 0, 0, 0, 0, 0, eigVal(4), 0, 0];
       [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, eigVal(5), 0];
       [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, eigVal(6)];];
   
% Reordering Eigen Value   
V_inv = [eigVec(:,7),eigVec(:,8),eigVec(:,9),eigVec(:,10),eigVec(:,11),eigVec(:,12),...
         eigVec(:,1),eigVec(:,2),eigVec(:,3),eigVec(:,4),eigVec(:,5),eigVec(:,6)];
%eigVal;

V11 = [[V_inv(1,1), V_inv(1,2), V_inv(1,3),V_inv(1,4), V_inv(1,5), V_inv(1,6)];
       [V_inv(2,1), V_inv(2,2), V_inv(2,3),V_inv(2,4), V_inv(2,5), V_inv(2,6)];
       [V_inv(3,1), V_inv(3,2), V_inv(3,3),V_inv(3,4), V_inv(3,5), V_inv(3,6)];
       [V_inv(4,1), V_inv(4,2), V_inv(4,3),V_inv(4,4), V_inv(4,5), V_inv(4,6)];
       [V_inv(5,1), V_inv(5,2), V_inv(5,3),V_inv(5,4), V_inv(5,5), V_inv(5,6)];
       [V_inv(6,1), V_inv(6,2), V_inv(6,3),V_inv(6,4), V_inv(6,5), V_inv(6,6)]];
   
V21 = [[V_inv(7,1), V_inv(7,2), V_inv(7,3),V_inv(7,4), V_inv(7,5), V_inv(7,6)];
       [V_inv(8,1), V_inv(8,2), V_inv(8,3),V_inv(8,4), V_inv(8,5), V_inv(8,6)];
       [V_inv(9,1), V_inv(9,2), V_inv(9,3),V_inv(9,4), V_inv(9,5), V_inv(9,6)];
       [V_inv(10,1), V_inv(10,2), V_inv(10,3),V_inv(10,4), V_inv(10,5), V_inv(10,6)];
       [V_inv(11,1), V_inv(11,2), V_inv(11,3),V_inv(11,4), V_inv(11,5), V_inv(11,6)];
       [V_inv(12,1), V_inv(12,2), V_inv(12,3),V_inv(12,4), V_inv(12,5), V_inv(12,6)]];
   
P = V21*inv(V11);
P = real(P);
inner = R + B_tld_y.'*P*B_hat;
F = (inv(inner)*B_tld_y.'*P*A_hat + inv(R)*W_y.');
alpha = - F;

%% Finally, Construct policy function
% Order [k, 1, z, tl, tx, g_ss]

% Constructing Policy function
k_grid = linspace(0.75*kss, 1.25*kss, 100);
k_pol = zeros(100, 1);
h_pol = zeros(100, 1);

for i = 1:100
    k_pol(i) = alpha(1,:)*[k_grid(i); 1; z_ss; tl_ss; tx_ss; g_ss];
    h_pol(i) = alpha(2,:)*[k_grid(i); 1; z_ss; tl_ss; tx_ss; g_ss];

end

figure(1);
tiledlayout(1,2)
nexttile
plot(k_grid, k_pol, 'LineWidth', 2, 'DisplayName','Capital Policy function');
title('Capital Policy Function (LQ)');
xlabel('level of Capital Today')
ylabel('level of Capital Tomorrow')

nexttile
plot(k_grid, h_pol, 'LineWidth', 2, 'DisplayName','Labour Policy function');
title('Labour Policy Function (LQ)');
xlabel('level of Capital Today')
ylabel('level of Hours')
saveas(figure(1), 'opt_pol_LQ.jpg')

toc;
