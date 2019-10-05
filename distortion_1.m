function [D , D_1 , D_2] = distortion_1(Pr , codebook_user_1 , codebook_user_2 , T , delta , f)
summation = 0 ;
summation_1 = 0 ;
summation_2 = 0 ;
parfor i_1 = 1 : 2
    for i_2 = 1 : 2
        for j_1 = 1 : 2
            for j_2 = 1 : 2
                
                i = (i_2 - 1) * 2 + i_1 ;
                j = (j_1 - 1) * 2 + j_2 ;
                
                u_1_index = find (T(: , 2) == i_1) ;
                u_1 = T(u_1_index , 1) ;
                
                u_2_index = find (T(: , 4) == i_2) ;
                for u_2_i = 1 : length(u_2_index)
                    u_2 = T(u_2_index (u_2_i) , 3) ;
                    
                    
                    hold_var = u_2 - codebook_user_2(j_1 , : , u_1_index) ;
                    hold_var = hold_var(:) ;
                    
                    summation = summation + Pr(i , j) * delta ...
                        * sum(f(u_1_index , u_2_index(u_2_i)) .* ((u_1 - codebook_user_1(j_2 , : , u_2_index(u_2_i))) .^ 2 + (hold_var) .^ 2));
                    
                    summation_1 = summation_1 + Pr (i , j) * delta...
                        * sum (f(u_1_index , u_2_index(u_2_i)) .* ((u_1 - codebook_user_1(j_2 , : , u_2_index(u_2_i))) .^ 2)) ;
                    summation_2 = summation_2 + Pr (i , j) * delta * sum(f(u_1_index , u_2_index(u_2_i)) .* ((hold_var) .^ 2)) ;
                end
            end
        end
    end
end
D = summation ;
D_1 = summation_1 ;
D_2 = summation_2 ;
end