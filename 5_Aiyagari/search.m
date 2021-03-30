function [low_adj, upp_adj, low_ind, upp_ind] = search(ass, edg, N_ass)
    srch_low = repmat(1, N_ass, 1);
    srch_high = repmat(50, N_ass, 1);
    
    for i = 1:N_ass
        if ass - edg(i) >= 0.0
                srch_low(i) = i;
        end
    end
    
    for j = 1:N_ass
        if ass - edg(j) < 0.0
                srch_high(j) = j;
        end
    end
    
    low_ind = max(srch_low);
    upp_ind = min(srch_high);
    
    low_adj = edg(low_ind);
    upp_adj = edg(upp_ind);
end

