
%% Iterating Riccati Difference Equation
function [P0, F0] = iter(P0, A_tld, B_tld, R, Q_tld, tol)
    iter = 0;
    A = R + B_tld.'*P0*B_tld;
    P1 = Q_tld + A_tld.'*P0*A_tld - A_tld.'*P0*B_tld*inv(A)*B_tld.'*P0*A_tld;
    Pnorm = norm(P0);
    pDiff = norm(P1-P0);
    
    F0 = inv(A)*B_tld.'*P0*A_tld;
    F1 = inv(A)*B_tld.'*P1*A_tld;
    Fnorm = norm(F0);
    fDiff = norm(F1-F0);
    
    while pDiff > tol*Pnorm && fDiff > tol*Fnorm
        A = R + B_tld.'*P1*B_tld;
        P2 = Q_tld + A_tld.'*P1*A_tld - A_tld.'*P1*B_tld*inv(A)*B_tld.'*P1*A_tld;
        Pnorm = norm(P1);
        pDiff = norm(P2-P1);
        
        F2 = inv(A)*B_tld.'*P2*A_tld;
        Fnorm = norm(F1);
        fDiff = norm(F2-F1);
        
        iter = iter + 1;
        P0 = P1;
        P1 = P2;
        F0 = F1;
        F1 = F2;
        
    end
    
        
end
