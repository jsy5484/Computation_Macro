function Y = find_Sb(X, ihat, hhat, yhat, ghat)
    % Inputs : logZt, tlt, txt, log(gt)
    
    
    global Params;
    
    z_ss = X(1); tl_ss = X(2); tx_ss = X(3); g_ss = log(mean(ghat));
    S_bar = [z_ss; tl_ss; tx_ss; g_ss];
    
    SS = BCA_findSS(S_bar);
    kss = SS(1); hss = SS(3); xss = SS(4); yss = SS(5);
    
    A = [kss, 0, hss, xss, yss];
    B = [mean(ihat)/(Params.grate - 1 + Params.delta), 0, mean(hhat), mean(ihat), mean(yhat)];

    Y = norm(A-B);
end