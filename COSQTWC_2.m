function [SDR , SDR_1 , SDR_2 , T , codebook_user_1 , codebook_user_2] = COSQTWC_2(f , p , T_u_1 , T_u_2 , Pr_z_1i_given_z_1i_1 , Pr_z_2i_given_z_2i_1 , P_bit_1 , width ,  codebook_user_1 , codebook_user_2)
delta = width (2) * width (1) ;

Threshold = 0.001 ;
D = [1 2] ;
FileID = fopen ('Results.txt', 'a') ;
while abs(D(2) - D(1)) / D(2) > (Threshold)
    D (1) = D(2) ;
    %% Optimal Centroids for user 1
    for u_2_index = 1 : length(T_u_2)
        u_2 = T_u_2(u_2_index , 1) ;
        x_2_1 = T_u_2(u_2_index , 2) ;
        for y_2_1 = 1 : 2
            for y_2_2 = 1 : 2
                y_2 = (y_2_1 - 1) * 2 + y_2_2 ;
                numerator = 0 ;
                denominator = 0 ;
                for y_1_1 = 1 : 2
                    x_2_2 = T_u_2(u_2_index , 2 + y_1_1) ;
                    
                     hold_var = (y_1_1 - 1 ) * 2 + y_2_1 ;
                     f_u_1_u_2_given_y_11_y_21 = f (: , : , hold_var) ;
                    for y_1_2 = 1 : 2
                        for x_1_2 = 1 : 2
                            
                            u_1_index = find (T_u_1(: , 2 + y_2_1) == x_1_2) ;
                            for u_1_i = 1 : length(u_1_index)
                                x_1_1 = T_u_1(u_1_index(u_1_i) , 2) ;
                                u_1 = T_u_1(u_1_index(u_1_i) , 1) ;
                                
                                binary_z_11 = xor(x_2_1 - 1 , y_1_1 - 1 ) ;
                                binary_z_12 = xor(x_2_2 - 1 , y_1_2 - 1) ;
                                
                                binary_z_21 = xor(x_1_1 - 1 , y_2_1 - 1) ;
                                binary_z_22 = xor(x_1_2 - 1 , y_2_2 - 1) ;
                                
                                
                                
                                numerator = numerator + P_bit_1(y_1_1 , y_2_1) ...
                                    * Pr_z_1i_given_z_1i_1(binary_z_11 + 1 , binary_z_12 + 1) ...
                                    * Pr_z_2i_given_z_2i_1(binary_z_21 + 1 , binary_z_22 + 1) ...
                                    * u_1 .* f_u_1_u_2_given_y_11_y_21(u_1_index(u_1_i) , u_2_index) ;
                                
                                denominator = denominator + P_bit_1(y_1_1 , y_2_1) ...
                                    * Pr_z_1i_given_z_1i_1(binary_z_11 + 1 , binary_z_12 + 1) ...
                                    * Pr_z_2i_given_z_2i_1(binary_z_21 + 1 , binary_z_22 + 1) ...
                                    * f_u_1_u_2_given_y_11_y_21(u_1_index(u_1_i) , u_2_index) ;
                            end
                        end
                    end
                end
                if (numerator == 0 && denominator == 0)
                    codebook_user_1(y_2, : , u_2_index) = p * u_2 ;
                else
                    codebook_user_1(y_2 , : , u_2_index) = numerator / denominator ;
                end
            end
        end
    end
    
    %% Optimal Centroids for user 2
    for u_1_index = 1 : length(T_u_1)
        u_1 = T_u_1(u_1_index , 1) ;
        x_1_1 = T_u_1(u_1_index , 2) ;
        for y_1_1 = 1 : 2
            for y_1_2 = 1 : 2
                y_1 = (y_1_1 - 1) * 2 + y_1_2 ;
                numerator = 0 ;
                denominator = 0 ;
                for y_2_1 = 1 : 2
                    x_1_2 = T_u_1(u_1_index , 2 + y_2_1) ;
                     hold_var = (y_1_1 - 1 ) * 2 + y_2_1 ;
                     f_u_1_u_2_given_y_11_y_21 = f (: , : , hold_var) ;
                    for y_2_2 = 1 : 2
                        for x_2_2 = 1 : 2
                            
                            
                            u_2_index = find(T_u_2(: , 2 + y_1_1) == x_2_2) ;
                            for u_2_i = 1 : length(u_2_index)
                                u_2 = T_u_2(u_2_index(u_2_i) , 1) ;
                                x_2_1 = T_u_2(u_2_index(u_2_i) , 2) ;
                                
                                binary_z_11 = xor(x_2_1 - 1 , y_1_1 - 1 ) ;
                                binary_z_12 = xor(x_2_2 - 1 , y_1_2 - 1) ;
                                
                                binary_z_21 = xor(x_1_1 - 1 , y_2_1 - 1) ;
                                binary_z_22 = xor(x_1_2 - 1 , y_2_2 - 1) ;
                                
                                
                                numerator = numerator + P_bit_1(y_1_1 , y_2_1) ...
                                    * Pr_z_1i_given_z_1i_1(binary_z_11 + 1 , binary_z_12 + 1) ...
                                    * Pr_z_2i_given_z_2i_1(binary_z_21 + 1 , binary_z_22 + 1) ...
                                    * u_2 .* f_u_1_u_2_given_y_11_y_21(u_1_index , u_2_index(u_2_i))' ;
                                denominator = denominator + P_bit_1(y_1_1 , y_2_1) ...
                                    * Pr_z_1i_given_z_1i_1(binary_z_11 + 1 , binary_z_12 + 1) ...
                                    * Pr_z_2i_given_z_2i_1(binary_z_21 + 1 , binary_z_22 + 1) ...
                                    * f_u_1_u_2_given_y_11_y_21(u_1_index , u_2_index(u_2_i))' ;
                            end
                        end
                    end
                end
                if (numerator == 0 && denominator == 0)
                    codebook_user_2(y_1 , : , u_1_index) = u_1 * p ;
                else
                    codebook_user_2(y_1 , : , u_1_index) = numerator / denominator ;
                end
            end
        end
    end
    
    %% Optimal Partitions for user 1
    parfor u_1_index = 1 : length(T_u_1)
        u_1 = T_u_1(u_1_index , 1) ;
        x_1_1 = T_u_1(u_1_index , 2) ;
        summation = 0 ;
        temp = zeros(2 , 2) ;
        for y_2_1 = 1 : 2
            for x_1_2 = 1 : 2
                for y_1_1 = 1 : 2
                     hold_var = (y_1_1 - 1 ) * 2 + y_2_1 ;
                     f_u_1_u_2_given_y_11_y_21 = f (: , : , hold_var) ;
                    for y_2_2 = 1 : 2
                        y_2 = (y_2_1 - 1) * 2 + y_2_2 ;
                        for y_1_2 = 1 : 2
                            y_1 = (y_1_1 - 1) * 2 + y_1_2 ;
                            for x_2_2 = 1 : 2
                                
                                u_2_index = find (T_u_2(: , 2 + y_1_1) == x_2_2) ;
                                
                                for u_2_i = 1 : length(u_2_index)
                                    x_2_1 = T_u_2(u_2_index(u_2_i) , 2) ;
                                    u_2 = T_u_2(u_2_index(u_2_i) , 1) ;
                                    
                                    binary_z_11 = xor(x_2_1 - 1 , y_1_1 - 1 ) ;
                                    binary_z_12 = xor(x_2_2 - 1 , y_1_2 - 1) ;
                                    
                                    binary_z_21 = xor(x_1_1 - 1 , y_2_1 - 1) ;
                                    binary_z_22 = xor(x_1_2 - 1 , y_2_2 - 1) ;
                                    
                                    
                                    hold_var = u_1 - codebook_user_1(y_2 , : , u_2_index(u_2_i)) ;
                                    hold_var = hold_var(:) ;
                                    
                                    summation = summation + P_bit_1(y_1_1 , y_2_1) ...
                                        * Pr_z_1i_given_z_1i_1(binary_z_11 + 1 , binary_z_12 + 1) ...
                                        * Pr_z_2i_given_z_2i_1(binary_z_21 + 1 , binary_z_22 + 1) ...
                                        * width(2) * f_u_1_u_2_given_y_11_y_21(u_1_index , u_2_index(u_2_i))' ...
                                        .* ((hold_var).^ 2 + (u_2 - codebook_user_2(y_1 , : , u_1_index)).^2 );
                                end
                            end
                        end
                    end
                end
                temp (y_2_1 , x_1_2) = summation ;
                summation = 0 ;
            end
        end
        [~ , partition_index] = min(temp , [] , 2) ;
        myT_u_1(u_1_index , :) = partition_index ;
    end
    T_u_1(: , 3 : 4) = myT_u_1 ;
    
    %% Optimal Partitions for user 2
    parfor u_2_index = 1 : length(T_u_2)
        u_2 = T_u_2(u_2_index , 1) ;
        x_2_1 = T_u_2(u_2_index , 2) ;
        summation = 0 ;
        temp = zeros (2 , 2) ;
        for y_1_1 = 1 : 2
            for x_2_2 = 1 : 2
                for y_2_1 = 1 : 2
                     hold_var = (y_1_1 - 1 ) * 2 + y_2_1 ;
                     f_u_1_u_2_given_y_11_y_21 = f (: , : , hold_var) ;
                    for y_1_2 = 1 : 2
                        y_1 = (y_1_1 - 1) * 2 + y_1_2 ;
                        for y_2_2 = 1 : 2
                            y_2 = (y_2_1 - 1) * 2 + y_2_2 ;
                            for x_1_2 = 1 : 2
                                
                                u_1_index = find (T_u_1(: , 2 + y_2_1) == x_1_2)  ;
                                for u_1_i = 1 : length(u_1_index)
                                    x_1_1 = T_u_1(u_1_index(u_1_i) , 2) ;
                                    u_1 = T_u_1(u_1_index(u_1_i) , 1) ;
                                    
                                    binary_z_11 = xor(x_2_1 - 1 , y_1_1 - 1 ) ;
                                    binary_z_12 = xor(x_2_2 - 1 , y_1_2 - 1) ;
                                    
                                    binary_z_21 = xor(x_1_1 - 1 , y_2_1 - 1) ;
                                    binary_z_22 = xor(x_1_2 - 1 , y_2_2 - 1) ;
                                    
                                    
                                    
                                    
                                    hold_var = u_2 - codebook_user_2(y_1 , : , u_1_index(u_1_i)) ;
                                    hold_var = hold_var(:) ;
                                    
                                    summation = summation + P_bit_1(y_1_1 , y_2_1) ...
                                        * Pr_z_1i_given_z_1i_1(binary_z_11 + 1 , binary_z_12 + 1) ...
                                        * Pr_z_2i_given_z_2i_1(binary_z_21 + 1 , binary_z_22 + 1) ...
                                        * width(1) * f_u_1_u_2_given_y_11_y_21(u_1_index(u_1_i) , u_2_index) ...
                                        .* ((u_1 - codebook_user_1(y_2 , : , u_2_index)).^ 2 + (hold_var).^2);
                                end
                            end
                        end
                    end
                end
                temp(y_1_1 , x_2_2) = summation ;
                summation = 0 ;
            end
        end
        [~ , partition_index] = min(temp , [] , 2) ;
        myT_u_2(u_2_index, :) = partition_index ;
    end
    T_u_2(: , 3 : 4) = myT_u_2 ;
    
    %% Distortion
    [D(2) , D_1 , D_2] = distortion_2(f , codebook_user_1 , codebook_user_2 , Pr_z_1i_given_z_1i_1 , Pr_z_2i_given_z_2i_1 , P_bit_1 , T_u_1 , T_u_2 , delta);
    fprintf (FileID , '\n Overall: \n') ;
    fprintf (FileID , ' %f ' , D(2)) ;
    fprintf (FileID , ' %f ' , D_1) ;
    fprintf (FileID , ' %f ' , D_2) ;
end
Distortion = D(2) ;
SDR = 10 * log10(2 / D(2)) ;
SDR_1 = 10 * log10(1 / D_1) ;
SDR_2 = 10 * log10(1 / D_2) ;
T = cat(2 , T_u_1 , T_u_2) ; 
fclose (FileID) ;
end