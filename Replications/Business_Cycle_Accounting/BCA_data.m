function Yhat = BCA_data(gamma_z)
temp = (1+gamma_z);
% Import Actual Data
Y = bringDat;
%importdata('data_fin.csv', ',');
prod = Y(:,1); invest = Y(:,2); hrs = Y(:,3); gov = Y(:,4);
T = length(prod);
%gz = 1.3791058031501657;
%gz = (1+gamma_z)^125;
% Apply Base Year Conversion
%bprod   =   prod/Y(125,1)*gz;
%binvest = invest/Y(125,1)*gz;
%bgov    =    gov/Y(125,1)*gz;

gz = (1+gamma_z)^81;
% Apply Base Year Conversion
bprod   =   prod/Y(81,1)*gz;
binvest = invest/Y(81,1)*gz;
bgov    =    gov/Y(81,1)*gz;


% Detrend the Data
yhat = zeros(T,1); ihat = zeros(T,1); hhat = zeros(T,1); ghat = zeros(T,1);
for i = 1:T
    yhat(i) = bprod(i)  /(temp^(i-1));
    ihat(i) = binvest(i)/(temp^(i-1));
    hhat(i) = hrs(i);
    ghat(i) = bgov(i)   /(temp^(i-1));
end

yhat = yhat(81:end); ihat = ihat(81:end); hhat = hhat(81:end); 
ghat = ghat(81:end);

%yhat = yhat(165:232); ihat = ihat(165:232); hhat = hhat(165:232); 
%ghat = ghat(165:232);

Yhat = [yhat, ihat, hhat, ghat];

end
