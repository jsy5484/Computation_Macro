
nsim = 500;         
nos = 74;         
T = 10000;        

for s = 1:nsim

    muz = 0;
    y_shock = rand(T,1);   
    zpath = zeros(T,1);
    zpath(1,1) = find(Z==max(Z(find(abs(muz - Z)==min(abs(muz - Z)))))); 

    for t = 2:T
        temp = cumsum(Zprob1(zpath(t - 1),:));
        if y_shock(t) < temp(1)
            zpath(t,1) = 1;
        elseif y_shock(t) > temp(ny)
            zpath(t,1) = ny;
        else
            zpath(t,1) = find(y_shock(t)<=temp(2:ny) & y_shock(t)>temp(1:ny-1))+1; 
        end
    end

    % Generate Redemption (1=redemption,0=autarky)
    redemp = ones(floor(theta*T),1);
    redemp = [redemp; zeros(T-length(redemp),1)];
    temp = randperm(T)'; redemp = redemp(temp);

    % Starting value
    bpath = zeros(T,1);
    defpath = zeros(T,1);  % Defaultpath(t)=1 if country defaulted in period t
    history = zeros(T,1);  % History(t)=1 if country is in autarky in period t.
    state = zeros(T,1);
    qpath = zeros(T,1);    % Vector of q 
    Edefpath = ones(T,1);  % Vector of default probabilities
    bpath(1) = bzero;      % Economies start without debt
    Vpath = zeros(T,1);    % Vector of value functions


    for t = 1:T-1
        state(t) = (zpath(t)-1)*nB + bpath(t);
        if history(t) == 0 || (history(t)==1 && redemp(t)==1) % if never defaulted or defaulted and redemption
            if default(state(t)) == 0                       % if no default in t
                bpath(t+1) = policygood(state(t));
                qpath(t) = q0(state(t),policygood(state(t)));
                Edefpath(t) = Edef(state(t),policygood(state(t)));
                Vpath(t) = Vc(state(t));
            elseif default(state(t)) == 1                   % if default in t
                bpath(t+1) = bzero;
                history(t+1) = 1;
                defpath(t) = 1;
                qpath(t) = 1/(1+r);
                Edefpath(t) = 0;
                Vpath(t) = Vd(state(t));
            end
            
        elseif history(t) == 1 && redemp(t) == 0 % if defaulted before
                bpath(t+1) = bzero;
                history(t+1) = 1;
                qpath(t) = 1/(1+r);
                Edefpath(t) = 0;
                Vpath(t) = Vd(state(t));
        end
    end

    qpath(T)=NaN;

    % Variables path
    ypath = exp(Z(zpath)); 
    b_y=b(bpath(2:T))./(ypath(1:T-1));

    bonds = b(bpath);
    dbonds = bonds(2:T)-bonds(1:T-1);
    dbondsy = dbonds./ypath(1:T-1);

    nx = bonds(2:T).*qpath(1:T-1)-bonds(1:T-1);
    nx(find(defpath(1:T-1)==1))=0;

    c = ypath(1:T-1) - nx;
    c(find(defpath(1:T-1)==1))=ypath(find(defpath(1:T-1)==1));

    sovr = 1./qpath - 1;
    spread = sovr-r;


    defaultt = find(defpath(1:T-1)==1); %find events of default
    defpath1=defpath;
    ssize=size(defaultt);


    % Guarantees that:(1) default event has at least #obs (=74) and
    %                 (2) default event has only 1 default in the period
    %                 (3) get only one random default event for each simulation
    for i=2:ssize
        if defaultt(1)<nos+1          %nos=74 (Arellano(2008))
           defpath1(defaultt(1))=0;
        end
        if defaultt(i)<nos+1
           defpath1(defaultt(i))=0;
        end
        if (defaultt(i)-defaultt(i-1))<nos+1
           defpath1(defaultt(i))=0;
        end       
    end

    % New default dummy for only the desired default events
    defaultt1 = find(defpath1(1:T-1)==1);
    ssize1=size(defaultt1);
    aleat=randi([1 ssize1(1)],1,1);

    %Variables
    %Output
    yt=ypath(defaultt1(aleat)-nos:defaultt1(aleat)-1);
    
    %Trade Balance
    nxt = nx(defaultt1(aleat)-nos:defaultt1(aleat)-1);
    nxy = nxt./yt;
    
    %Consumption
    ct=c(defaultt1(aleat)-nos:defaultt1(aleat)-1);
    
    %Spread
    spreadt=spread(defaultt1(aleat)-nos:defaultt1(aleat)-1);
    spreadtannual=(1+spreadt).^4-1;
    
    %Bonds
    Bt=bonds(defaultt1(aleat)-nos:defaultt1(aleat)-1);
    By=Bt./yt;
    dbondsyt=dbondsy(defaultt1(aleat)-nos:defaultt1(aleat)-1);
    
    %Correlations
    STD=std([yt, spreadtannual, nxy, ct]);
    CC=corrcoef([yt, spreadtannual, nxy, ct]);
    AC=corrcoef(yt(1:end-1),yt(2:end));
    stdyrnxc(s,:)=STD;
    ccr(s,:)=[CC(1,2:4),CC(2,3:4)];
    defaultpc(s)=mean(defpath(find(bpath<bzero)))*100;
    
    y_sim(s,:)=yt;
    spread_sim(s,:)=spreadtannual;
    nxy_sim(s,:)=nxy;
    c_sim(s,:)=ct;
    by_sim(s,:)=By;
  
end

%Mean statistics
meansp=mean(spread_sim)*100;
meannxy=mean(nxy_sim)*10;
meanc=mean(c_sim)*10;
meany=mean(y_sim)*10;
meanBy=mean(by_sim)*100;

disp('Default_event Spread Trade_Balance Consumption Output Debt/Output')
disp([meansp(end),meannxy(end),meanc(end),meany(end), meanBy(end)])

disp('mean std Output Spread Trade_Balance Consumption')
disp(mean(stdyrnxc,1))

disp('mean correlation Output-Spread Output-Trade_Balance Output-Consumption Spread-Trade_Balance Output-Output')
disp(mean(ccr,1))
disp('mean default')
disp(mean(defaultpc))

figure6=figure(6);
titulo=sprintf('Default events');
plot(defpath)
xlabel('Time')
ylabel('Default events')
title(titulo)


figure7=figure(7); 
[hAx,hLine1,hLine2] = plotyy(1:nos,meany(1:nos),1:nos-3,meannxy(1:nos-3));

title('Output x Trade Balance')
xlabel('Time')

ylabel(hAx(1),'Output') 
ylabel(hAx(2),'Trade Balance') 

figure8=figure(8); 
[hAx,hLine1,hLine2] = plotyy(1:nos,meany,1:nos,meansp);

title('Output x Spread')
xlabel('Time')

ylabel(hAx(1),'Output') 
ylabel(hAx(2),'Spread') 

figure9=figure(9);
titulo=sprintf('Bonds over Output');
plot(meanBy)
xlabel('Time')
ylabel('Bonds over Y')
title(titulo)

figure10=figure(10);
titulo=sprintf('Output');
plot(meany)
xlabel('Time')
ylabel('meany')
title(titulo)

figure11=figure(11);
titulo=sprintf('Trade Balance over Output');
plot(1:nos-3,meannxy(1:nos-3))
xlabel('Time')
ylabel('meannxy')
title(titulo)

figure12=figure(12);
titulo=sprintf('Consumption');
plot(meanc)
xlabel('Time')
ylabel('meanc')
title(titulo)

figure13=figure(13);
titulo=sprintf('Spread');
plot(meansp)
xlabel('Time')
ylabel('meansp')
title(titulo)



