function [KM_time, KM_cross]  = aggSimul(horiz, idio_shock, agg_shock, KM_max, KM_min, k_pol, KM_grid, k_grid, k_min, k_max,KM_cross)
% Using policy function, simulate cross-sectional distribution, and
% record first-moments and distribution every periods.

KM_time = zeros(horiz,1); a2 = [1; 2]; epsilon2 = [1; 2];

    for t = 1:horiz
       KM_time(t) = mean(KM_cross); 
       KM_time(t) = KM_time(t)*(KM_time(t)>=KM_min)*(KM_time(t)<=KM_max)+ ...
                    KM_min*(KM_time(t)<KM_min)+KM_max*(KM_time(t)>KM_max); 

       k_pol_dim = interpn(k_grid,    KM_grid, a2,           epsilon2, k_pol, ...
                           k_grid, KM_time(t), agg_shock(t), epsilon2,'spline');
       k_pol_t = squeeze(k_pol_dim); 


       KM_cross_new = interpn(k_grid, epsilon2, k_pol_t, ...
                              KM_cross, idio_shock(t,:),'spline'); 

       KM_cross_new = KM_cross_new.*(KM_cross_new >= k_min).*(KM_cross_new<=k_max)+ ...
                     k_min*(KM_cross_new<k_min)+k_max*(KM_cross_new>k_max);               
       KM_cross = KM_cross_new;

    end            
   
end