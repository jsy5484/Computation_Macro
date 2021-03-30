function [agg_shock, idy_shock] = shock_simul(trans, NH, horiz, unemp_b)
    agg_shock = zeros(horiz,1);   % vector of aggregate shocks
    idy_shock = zeros(horiz,NH);  % matrix of idiosyncratic shocks 
    
    prob_agg = zeros(2,2);  
    
    prob_agg(1,1) = trans(1,1) + trans(1,2); 
    prob_agg(2,1) = 1 - prob_agg(1,1);  
    prob_agg(2,2) = trans(3,3) + trans(3,4); 
    prob_agg(1,2) = 1 - prob_agg(2,2);

    pr_bb_uu = trans(1,1)/prob_agg(1,1); 
    pr_bb_ue = 1 - pr_bb_uu;
    
    pr_bb_ee = trans(2,2)/prob_agg(1,1); 
    pr_bb_eu = 1 - pr_bb_ee;
    
    pr_bg_uu = trans(1,3)/prob_agg(2,1); 
    pr_bg_ue = 1 - pr_bg_uu;
    
    pr_bg_ee = trans(2,4)/prob_agg(2,1); 
    pr_bg_eu = 1 - pr_bg_ee;
    
    pr_gb_uu = trans(3,1)/prob_agg(1,2); 
    pr_gb_ue = 1 - pr_gb_uu;
    
    pr_gb_ee = trans(4,2)/prob_agg(1,2); 
    pr_gb_eu = 1 - pr_gb_ee;
    
    pr_gg_uu = trans(3,3)/prob_agg(2,2); 
    pr_gg_ue = 1 - pr_gg_uu;
    
    pr_gg_ee = trans(4,4)/prob_agg(2,2); 
    pp_gg_eu = 1 - pr_gg_ee;
    
    
    %% PART 1. Construct Simulation of Aggregate Shock
    agg_shock(1) = 1; % assume initial aggregate state is bad 
    
    for t = 2:horiz
       z = rand; 
       if z <= prob_agg(1,agg_shock(t-1)) 
          agg_shock(t) = 1; 
       else
          agg_shock(t) = 2;
       end
    end
    
    %% PART 2. Construct Simulation of Idiosyncratic Shock
    
    % Draw an Initial Distribution
    for i = 1:NH
        z = rand;
        if z <= unemp_b 
            idy_shock(1,i)=1;
        else
            idy_shock(1,i)=2;
        end
    end

    % Afterwards
    for t = 2:horiz
        % Case 1  (Bad to Bad)
        if agg_shock(t-1) == 1 && agg_shock(t) == 1 
            for i = 1:NH
                z = rand;
                if idy_shock(t-1,i)==1 
                    if z <= pr_bb_uu
                       idy_shock(t,i)=1;
                    else
                       idy_shock(t,i)=2;
                    end
                 else                
                     if z <= pr_bb_ee
                       idy_shock(t,i)=2;
                    else
                       idy_shock(t,i)=1;
                     end
                end
            end
        end
       % Case 2  (Bad to Good)
       if agg_shock(t-1)==1 && agg_shock(t)==2 
          for i = 1:NH
             z = rand;
             if idy_shock(t-1,i)==1
                if z <= pr_bg_uu
                   idy_shock(t,i)=1;
                else
                   idy_shock(t,i)=2;
                end
             else
                if z <= pr_bg_ee
                   idy_shock(t,i)=2;
                else
                   idy_shock(t,i)=1;
                end
             end
          end
       end
       % Case 3  (Good to Bad)
       if agg_shock(t-1)==2 && agg_shock(t)==1 
          for i = 1:NH
             z = rand;
             if idy_shock(t-1,i)==1
                if z <= pr_gb_uu
                   idy_shock(t,i)=1;
                else
                   idy_shock(t,i)=2;
                end
             else
                if z <= pr_gb_ee
                   idy_shock(t,i)=2;
                else
                   idy_shock(t,i)=1;
                end
             end
          end
       end
       % Case 4 (Good to Good)
       if agg_shock(t-1)==2 && agg_shock(t)==2 
          for i=1:NH
             z =rand;
             if idy_shock(t-1,i)==1
                if z <= pr_gg_uu
                   idy_shock(t,i)=1;
                else
                   idy_shock(t,i)=2;
                end
             else
                if z <= pr_gg_ee
                   idy_shock(t,i)=2;
                else
                   idy_shock(t,i)=1;
                end
             end
          end
       end
    
    end
end

