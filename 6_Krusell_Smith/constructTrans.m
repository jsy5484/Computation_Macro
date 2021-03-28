function transition = constructTrans(dur_ug, dur_ub, dur_gd, dur_bd, unemp_gd, unemp_bd)
% Construct a transition probability from information of aggregate economy 

    % For aggregate Shock
    p_gg = (dur_gd-1)/dur_gd;
    p_gb = 1 - (dur_bd-1)/dur_bd;
    p_bg = 1 - (dur_gd-1)/dur_gd;
    p_bb = (dur_bd-1)/dur_bd;

    p_gg00 = (dur_ug-1)/dur_ug;
    p_bb00 = (dur_ub-1)/dur_ub;
    p_bg00 = 1.25*p_bb00;
    p_gb00 = 0.75*p_gg00;
    p_gg01 = (unemp_gd - unemp_gd*p_gg00)/(1-unemp_gd);
    p_bb01 = (unemp_bd - unemp_bd*p_bb00)/(1-unemp_bd);
    p_bg01 = (unemp_bd - unemp_gd*p_bg00)/(1-unemp_gd);
    p_gb01 = (unemp_gd - unemp_bd*p_gb00)/(1-unemp_bd);

    p_gg10 = 1 - (dur_ug-1)/dur_ug;
    p_bb10 = 1 - (dur_ub-1)/dur_ub;
    p_bg10 = 1 - 1.25*p_bb00;
    p_gb10 = 1 - 0.75*p_gg00;
    p_gg11 = 1 - (unemp_gd - unemp_gd*p_gg00)/(1-unemp_gd);
    p_bb11 = 1 - (unemp_bd - unemp_bd*p_bb00)/(1-unemp_bd);
    p_bg11 = 1 - (unemp_bd - unemp_gd*p_bg00)/(1-unemp_gd);
    p_gb11 = 1 - (unemp_gd - unemp_bd*p_gb00)/(1-unemp_bd);

    % Joint transition matrix
    Pr_11 = p_gg*p_gg11;
    Pr_21 = p_bg*p_bg11;
    Pr_31 = p_gg*p_gg01;
    Pr_41 = p_bg*p_bg01;

    Pr_12 = p_gb*p_gb11;
    Pr_22 = p_bb*p_bb11;
    Pr_32 = p_gb*p_gb01;
    Pr_42 = p_bb*p_bb01;

    Pr_13 = p_gg*p_gg10;
    Pr_23 = p_bg*p_bg10;
    Pr_33 = p_gg*p_gg00;
    Pr_43 = p_bg*p_bg00;

    Pr_14 = p_gb*p_gb10;
    Pr_24 = p_bb*p_bb10;
    Pr_34 = p_gb*p_gb00;
    Pr_44 = p_bb*p_bb00;

    transition =[[Pr_11, Pr_21, Pr_31, Pr_41];
                 [Pr_12, Pr_22, Pr_32, Pr_42];
                 [Pr_13, Pr_23, Pr_33, Pr_43];
                 [Pr_14, Pr_24, Pr_34, Pr_44]];

end