function [SDR , SDR_1 , SDR_2 , T , codebook_user_1 , codebook_user_2] = COSQTWC_4(p , f , T_u_1 , T_u_2 , Pr_z_1i_given_z_1i_1 , Pr_z_2i_given_z_2i_1 , P_bit_3 , width ,  codebook_user_1 , codebook_user_2)
delta = width(2) * width(1);


Threshold = 0.001 ;
D = [1 2] ;
FileID = fopen ('Results.txt', 'a') ;
while abs(D(2) - D(1)) / D(2) > (Threshold)
    D (1) = D(2) ;
    
    %% Optimal Centroids for user 1
    for u_2_index = 1 : length(T_u_2)
        u_2 = T_u_2(u_2_index , 1) ;
        for y_2_1 = 1 : 2
            for y_2_2 = 1 : 2
                y_2_12 = (y_2_1 - 1) * 2 + y_2_2 ;
                for y_2_3 = 1 : 2
                    y_2_prime = (y_2_1 - 1) * 4 + (y_2_2 - 1) * 2 + y_2_3 ;
                    for y_2_4 = 1 : 2
                        y_2 = (y_2_1 - 1) * 8 + (y_2_2 - 1) * 4 + (y_2_3 - 1) * 2 + y_2_4 ;
                        numerator = 0 ;
                        denominator = 0 ;
                        for y_1_1 = 1 : 2
                            for y_1_2 = 1 : 2
                                y_1_12 = (y_1_1 - 1) * 2 + y_1_2 ;
                                x_2_3 = T_u_2(u_2_index , 1 + y_1_12) ;
                                for y_1_3 = 1 : 2
                                    y_1_prime = (y_1_1 - 1) * 4 + (y_1_2 - 1) * 2 + y_1_3 ;
                                    x_2_4 = T_u_2(u_2_index , 5 + y_1_prime) ;
                                    for x_1_4 = 1 : 2
                                        for y_1_4 = 1 : 2
                                            
                                            hold_var = (y_1_1 - 1 ) * 32 + (y_1_2 - 1) * 16 + (y_1_3 - 1) * 8 + (y_2_1 - 1) * 4 + (y_2_2 - 1) * 2 + y_2_3;
                                            f_u_1_u_2_given_y_11_y_12_y_13_y_21_y_22_y_23 = f(: , : , hold_var) ;
                                            
                                            u_1_index = find (T_u_1(:, 5 + y_2_prime) == x_1_4) ;
                                            for u_1_i = 1 : length(u_1_index)
                                                u_1 = T_u_1(u_1_index(u_1_i) , 1) ;
                                                x_1_3 = T_u_1(u_1_index(u_1_i) , 1 + y_2_12) ;
                                                
                                                binary_z_13 = xor(x_2_3 - 1 , y_1_3 - 1 ) ;
                                                binary_z_14 = xor(x_2_4 - 1 , y_1_4 - 1) ;
                                                
                                                binary_z_23 = xor(x_1_3 - 1 , y_2_3 - 1) ;
                                                binary_z_24 = xor(x_1_4 - 1 , y_2_4 - 1) ;
                                                
                                                
                                                numerator = numerator + P_bit_3(y_1_prime , y_2_prime) ...
                                                    * Pr_z_1i_given_z_1i_1(binary_z_13 + 1 , binary_z_14 + 1) ...
                                                    * Pr_z_2i_given_z_2i_1(binary_z_23 + 1 , binary_z_24 + 1) ...
                                                    * u_1 .* f_u_1_u_2_given_y_11_y_12_y_13_y_21_y_22_y_23(u_1_index(u_1_i) , u_2_index) ;
                                                
                                                denominator = denominator + P_bit_3(y_1_prime , y_2_prime) ...
                                                    * Pr_z_1i_given_z_1i_1(binary_z_13 + 1 , binary_z_14 + 1) ...
                                                    * Pr_z_2i_given_z_2i_1(binary_z_23 + 1 , binary_z_24 + 1) ...
                                                    * f_u_1_u_2_given_y_11_y_12_y_13_y_21_y_22_y_23(u_1_index(u_1_i) , u_2_index) ;
                                            end
                                        end
                                    end
                                end
                            end
                        end
                        if (numerator == 0 && denominator == 0)
                            codebook_user_1(y_2 , : , u_2_index) = p * u_2 ;
                        else
                            codebook_user_1(y_2 , : , u_2_index) = numerator / denominator ;
                        end
                    end
                end
            end
        end
    end
    %% Optimal Centroids for user 2
    for u_1_index = 1 : length(T_u_1)
        u_1 = T_u_1(u_1_index , 1) ;
        for y_1_1 = 1 : 2
            for y_1_2 = 1 : 2
                y_1_12 = (y_1_1 - 1) * 2 + y_1_2 ;
                
                for y_1_3 = 1 : 2
                    y_1_prime = (y_1_1 - 1) * 4 + (y_1_2 - 1) * 2 + y_1_3 ;
                    for y_1_4 = 1 : 2
                        y_1 = (y_1_1 - 1) * 8 + (y_1_2 - 1) * 4 + (y_1_3 - 1) * 2 + y_1_4 ;
                        numerator = 0 ;
                        denominator = 0 ;
                        for y_2_1 = 1 : 2
                            for y_2_2 = 1 : 2
                                y_2_12 = (y_2_1 - 1) * 2 + y_2_2 ;
                                x_1_3 = T_u_1(u_1_index , 1 + y_2_12) ;
                                for y_2_3 = 1 : 2
                                    y_2_prime = (y_2_1 - 1) * 4 + (y_2_2 - 1) * 2 + y_2_3 ;
                                    x_1_4 = T_u_1(u_1_index , 5 + y_2_prime) ;
                                    for x_2_4 = 1 : 2
                                        for y_2_4 = 1 : 2
                                            
                                            hold_var = (y_1_1 - 1 ) * 32 + (y_1_2 - 1) * 16 + (y_1_3 - 1) * 8 + (y_2_1 - 1) * 4 + (y_2_2 - 1) * 2 + y_2_3;
                                            f_u_1_u_2_given_y_11_y_12_y_13_y_21_y_22_y_23 = f(: , : , hold_var) ;
                                            
                                            u_2_index = find (T_u_2(: , 5 + y_1_prime) == x_2_4) ;
                                            for u_2_i = 1 : length(u_2_index)
                                                u_2 = T_u_2(u_2_index(u_2_i) , 1) ;
                                                x_2_3 = T_u_2(u_2_index(u_2_i) , 1 + y_1_12) ;
                                                binary_z_13 = xor(x_2_3 - 1 , y_1_3 - 1 ) ;
                                                binary_z_14 = xor(x_2_4 - 1 , y_1_4 - 1) ;
                                                
                                                binary_z_23 = xor(x_1_3 - 1 , y_2_3 - 1) ;
                                                binary_z_24 = xor(x_1_4 - 1 , y_2_4 - 1) ;
                                                
                                                
                                                
                                                numerator = numerator + P_bit_3(y_1_prime , y_2_prime)...
                                                    * Pr_z_1i_given_z_1i_1(binary_z_13 + 1 , binary_z_14 + 1) ...
                                                    * Pr_z_2i_given_z_2i_1(binary_z_23 + 1 , binary_z_24 + 1) ...
                                                    * u_2 .* f_u_1_u_2_given_y_11_y_12_y_13_y_21_y_22_y_23(u_1_index , u_2_index(u_2_i))' ;
                                                
                                                denominator = denominator + P_bit_3(y_1_prime , y_2_prime) ...
                                                    * Pr_z_1i_given_z_1i_1(binary_z_13 + 1 , binary_z_14 + 1) ...
                                                    * Pr_z_2i_given_z_2i_1(binary_z_23 + 1 , binary_z_24 + 1) ...
                                                    * f_u_1_u_2_given_y_11_y_12_y_13_y_21_y_22_y_23(u_1_index , u_2_index(u_2_i))' ;
                                            end
                                        end
                                    end
                                end
                            end
                        end
                        if (numerator == 0 && denominator == 0)
                            codebook_user_2 (y_1 , : , u_1_index ) = p * u_1 ;
                        else
                            codebook_user_2 (y_1 , : , u_1_index) = numerator / denominator ;
                        end
                    end
                end
            end
        end
    end
    
    %% Optimal Partitions for user 1
    parfor u_1_index = 1 : length(T_u_1)
        summation = 0 ;
        temp = zeros (8 , 2) ;
        u_1 = T_u_1(u_1_index ,  1) ;
        for y_2_1 = 1 : 2
            for y_2_2 = 1 : 2
                y_2_12 = (y_2_1 - 1) * 2 + y_2_2 ;
                x_1_3 = T_u_1(u_1_index , 1 + y_2_12) ;
                for y_2_3 = 1 : 2
                    y_2_prime = (y_2_1 - 1) * 4 + (y_2_2 - 1) * 2 + y_2_3 ;
                    for x_1_4 = 1 : 2
                        for y_1_1 = 1 : 2
                            for y_1_2 = 1 : 2
                                y_1_12 = (y_1_1 - 1) * 2 + y_1_2 ;
                                for y_1_3 = 1 : 2
                                    y_1_prime = (y_1_1 - 1) * 4 + (y_1_2 - 1) * 2 + y_1_3 ;
                                    for y_2_4 = 1 : 2
                                        y_2 = (y_2_1 - 1) * 8 + (y_2_2 - 1) * 4 + (y_2_3 - 1) * 2 + y_2_4 ;
                                        for y_1_4 = 1 : 2
                                            y_1 = (y_1_1 - 1) * 8 + (y_1_2 - 1) * 4 + (y_1_3 - 1) * 2 + y_1_4 ;
                                            for x_2_4 = 1 : 2
                                                
                                                hold_var = (y_1_1 - 1 ) * 32 + (y_1_2 - 1) * 16 + (y_1_3 - 1) * 8 + (y_2_1 - 1) * 4 + (y_2_2 - 1) * 2 + y_2_3;
                                                f_u_1_u_2_given_y_11_y_12_y_13_y_21_y_22_y_23 = f(: , : , hold_var) ;
                                                
                                                
                                                u_2_index = find (T_u_2(: , 5 + y_1_prime) == x_2_4) ;
                                                for u_2_i = 1 : length(u_2_index)
                                                    u_2 = T_u_2(u_2_index(u_2_i) , 1) ;
                                                    x_2_3 = T_u_2(u_2_index(u_2_i) , 1 + y_1_12) ;
                                                    
                                                    binary_z_13 = xor(x_2_3 - 1 , y_1_3 - 1 ) ;
                                                    binary_z_14 = xor(x_2_4 - 1 , y_1_4 - 1) ;
                                                    
                                                    binary_z_23 = xor(x_1_3 - 1 , y_2_3 - 1) ;
                                                    binary_z_24 = xor(x_1_4 - 1 , y_2_4 - 1) ;
                                                    
                                                    
                                                    
                                                    hold_var = u_1 - codebook_user_1(y_2 , : , u_2_index(u_2_i)) ;
                                                    hold_var = hold_var(:) ;
                                                    
                                                    summation = summation + P_bit_3(y_1_prime , y_2_prime) ...
                                                        * Pr_z_1i_given_z_1i_1(binary_z_13 + 1 , binary_z_14 + 1) ...
                                                        * Pr_z_2i_given_z_2i_1(binary_z_23 + 1 , binary_z_24 + 1) ...
                                                        * delta * f_u_1_u_2_given_y_11_y_12_y_13_y_21_y_22_y_23(u_1_index , u_2_index(u_2_i))' .* ((hold_var).^2 ...
                                                        + (u_2 - codebook_user_2(y_1 , : , u_1_index)).^2) ;
                                                end
                                            end
                                        end
                                    end
                                end
                            end
                        end
                        temp(y_2_prime , x_1_4) = summation ;
                        summation = 0 ;
                    end
                end
            end
        end
        [~ , partition_index] = min(temp , [] , 2) ;
        myT_u_1 (u_1_index , :) = partition_index ;
    end
    T_u_1(: , 6 : 13) = myT_u_1 ;
    
    %% Optimal Partitions for user 2
    parfor u_2_index = 1 : length(T_u_2)
        summation = 0 ;
        temp = zeros (8 , 2) ;
        u_2 = T_u_2(u_2_index , 1) ;
        for y_1_1 = 1 : 2
            for y_1_2 = 1 : 2
                y_1_12 = (y_1_1 - 1) * 2 + y_1_2 ;
                x_2_3 = T_u_2(u_2_index , 1 + y_1_12) ;
                for y_1_3 = 1 : 2
                    y_1_prime = (y_1_1 - 1) * 4 + (y_1_2 - 1) * 2 + y_1_3 ;
                    for x_2_4 = 1 : 2
                        for y_2_1 = 1 : 2
                            for y_2_2 = 1 : 2
                                y_2_12 = (y_2_1 - 1) * 2 + y_2_2 ;
                                for y_2_3 = 1 : 2
                                    y_2_prime = (y_2_1 - 1) * 4 + (y_2_2 - 1) * 2 + y_2_3 ;
                                    for y_1_4 = 1 : 2
                                        y_1 = (y_1_1 - 1) * 8 + (y_1_2 - 1) * 4 + (y_1_3 - 1) * 2 + y_1_4 ;
                                        for y_2_4 = 1 : 2
                                            y_2 = (y_2_1 - 1) * 8 + (y_2_2 - 1) * 4 +(y_2_3 - 1) * 2 + y_2_4 ;
                                            for x_1_4 = 1 : 2
                                                
                                                hold_var = (y_1_1 - 1 ) * 32 + (y_1_2 - 1) * 16 + (y_1_3 - 1) * 8 + (y_2_1 - 1) * 4 + (y_2_2 - 1) * 2 + y_2_3;
                                                f_u_1_u_2_given_y_11_y_12_y_13_y_21_y_22_y_23 = f(: , : , hold_var) ;
                                                
                                                u_1_index = find (T_u_1(: , 5 + y_2_prime) == x_1_4 ) ;
                                                for u_1_i = 1 : length(u_1_index)
                                                    u_1 = T_u_1(u_1_index(u_1_i) , 1) ;
                                                    x_1_3 = T_u_1(u_1_index(u_1_i) , 1 + y_2_12) ;
                                                    
                                                    binary_z_13 = xor(x_2_3 - 1 , y_1_3 - 1 ) ;
                                                    binary_z_14 = xor(x_2_4 - 1 , y_1_4 - 1) ;
                                                    
                                                    binary_z_23 = xor(x_1_3 - 1 , y_2_3 - 1) ;
                                                    binary_z_24 = xor(x_1_4 - 1 , y_2_4 - 1) ;
                                                    
                                                    
                                                    hold_var = u_2 - codebook_user_2(y_1 , : , u_1_index(u_1_i)) ;
                                                    hold_var = hold_var(:) ;
                                                    summation = summation + P_bit_3(y_1_prime , y_2_prime ) ...
                                                        * Pr_z_1i_given_z_1i_1(binary_z_13 + 1 , binary_z_14 + 1) ...
                                                        * Pr_z_2i_given_z_2i_1(binary_z_23 + 1 , binary_z_24 + 1) ...
                                                        * delta * f_u_1_u_2_given_y_11_y_12_y_13_y_21_y_22_y_23(u_1_index(u_1_i) , u_2_index) ...
                                                        .* ((u_1 - codebook_user_1(y_2 , : , u_2_index)) .^ 2 ...
                                                        + (hold_var).^2 ) ;
                                                end
                                                
                                            end
                                        end
                                    end
                                end
                            end
                        end
                        temp (y_1_prime , x_2_4) = summation ;
                        summation = 0 ;
                    end
                end
            end
        end
        [~ , partition_index] = min(temp , [] , 2) ;
        myT_u_2(u_2_index , : ) = partition_index ;
    end
    T_u_2 (: , 6 : 13) = myT_u_2 ;
    
    
    %% Distortion
    [D(2) , D_1 , D_2] = distortion_4 (f , P_bit_3 , ...
        Pr_z_1i_given_z_1i_1 , Pr_z_2i_given_z_2i_1 , ...
        T_u_1 , T_u_2 , codebook_user_1 , codebook_user_2 , delta) ;
    fprintf (FileID , '\n Overall: \n') ;
    fprintf (FileID , ' %f ' , D(2)) ;
    fprintf (FileID , ' %f ' , D_1) ;
    fprintf (FileID , ' %f ' , D_2) ;
end
SDR = 10 * log10(2 / D(2)) ;
SDR_1 = 10 * log10(1 / D_1) ;
SDR_2 = 10 * log10(1 / D_2) ;
T = cat(2 , T_u_1 , T_u_2) ; 
fclose (FileID) ;
end