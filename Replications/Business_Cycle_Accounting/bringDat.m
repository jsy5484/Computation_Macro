function Y = bringDat

data = readtable('data.csv');
rGDP = table2array(data(:,2));
rPCE = table2array(data(:,3));
pPCE = table2array(data(:,4));
pCD = table2array(data(:,5));
pCND = table2array(data(:,6));
pCS = table2array(data(:,7));
nCD = table2array(data(:,8));
nCND = table2array(data(:,9));
nCS = table2array(data(:,10));
rGPDI = table2array(data(:,11));
rEX = table2array(data(:,12));
rIM = table2array(data(:,13));
rG = table2array(data(:,14));
pGC = table2array(data(:,15));
nGC = table2array(data(:,16));
pGI = table2array(data(:,17));
nGI = table2array(data(:,18));
nFETX = table2array(data(:,19));
nSTX = table2array(data(:,20));
nKCD = table2array(data(:,21));
nDCD = table2array(data(:,22));
hours = table2array(data(:,23));
pop = table2array(data(:,24));
%rGC = table2array(data(:,25));

rCD = 100*nCD./pCD;
rCND = 100*nCND./pCND;
rCS = 100*nCS./pCS;
rGC = 100*nGC./pGC;
rGI = 100*nGI./pGI;
nTX = nFETX + nSTX;
rTX = 100*nTX./pPCE;
rKCD = 100*nKCD./pCD;
rDCD = 100*nDCD./pCD;

Y = rGDP-rTX+0.04*rKCD+rDCD;
C = rCND+rCS-(rCND+rCS)./(rCND+rCS+rCD).*rTX+0.04*rKCD+rDCD;
X = rCD+rGPDI+rGI-rCD./(rCND+rCS+rCD).*rTX;
G = rGC+rEX-rIM;

prd   = Y./hours*10^9;
hpc   = hours./pop;
ypc   = Y./pop*10^9;
cpc   = C./pop*10^9;
xpc   = X./pop*10^9;
gpc   = G./pop*10^9;

%gz = ((1.016)^(1/4))^81;
%ZVAR = [ypc/ypc(81)*gz, xpc/ypc(81)*gz, hpc/1300, gpc/ypc(81)*gz];
%ZVAR = [ypc, xpc, hpc./1300, gpc];
%T = length(ZVAR(1));
%for t = 1:T
%    Y = [log(ZVAR(:,1))-log((1+gamma_z)^(t-1)), ...
%         log(ZVAR(:,2)) - log((1+gamma_z)^(t-1)),...
%         log(ZVAR(:,3)), ...
%         log(ZVAR(:,4)) - log((1+gamma_z)^t)];
%end
%lyt      = ZVAR(:,1);
%lxt      = ZVAR(:,2);
%llt      = ZVAR(:,3);
%lgt      = ZVAR(:,4);

Y = [ypc, xpc, hpc./1300, gpc];

%El_dat_y = table2array(uszvarq1(:,2));
%El_dat_x = table2array(uszvarq1(:,3));
%El_dat_h = table2array(uszvarq1(:,4));
%El_dat_g = table2array(uszvarq1(:,5));


%figure(1)
%plot(linspace(1,236,236), lyt)
%hold on;
%plot(linspace(1,183,183), El_dat_y)
%hold off;

%figure(2)
%plot(linspace(1,236,236), lxt)
%hold on;
%plot(linspace(1,183,183), El_dat_x)
%hold off;

%figure(3)
%plot(linspace(1,236,236), llt)
%hold on;
%plot(linspace(1,183,183), El_dat_h)
%hold off;


%figure(4)
%plot(linspace(1,236,236), lgt)
%hold on;
%plot(linspace(1,183,183), El_dat_g)
%hold off;
