function [D , D_1 , D_2] = distortion_4(f , P_bit_3 , Pr_z_1i_given_z_1i_1 , Pr_z_2i_given_z_2i_1 , T_u_1 , T_u_2 , codebook_user_1 , codebook_user_2 , delta)
summation = 0 ;
summation_1 = 0 ;
summation_2 = 0 ;
parfor y_1_1 = 1 : 2
    for y_1_2 = 1 : 2
        y_1_12 = (y_1_1 - 1) * 2 + y_1_2 ;
        for y_1_3 = 1 : 2
            y_1_prime = (y_1_1 - 1) * 4 + (y_1_2 - 1) * 2 + y_1_3 ;
            for y_2_1 = 1 : 2
                for y_2_2 = 1 : 2
                    y_2_12 = (y_2_1 - 1) * 2 + y_2_2 ;
                    for y_2_3 = 1 : 2
                        y_2_prime = (y_2_1 - 1) * 4 + (y_2_2 - 1) * 2 + y_2_3 ;
                        for x_1_4 = 1 : 2
                            for x_2_4 = 1 : 2
                                for y_1_4 = 1 : 2
                                    y_1 = (y_1_1 - 1) * 8 + (y_1_2 - 1) * 4 + (y_1_3 - 1) * 2 + y_1_4 ;
                                    for y_2_4 = 1 : 2
                                        y_2 = (y_2_1 - 1) * 8 + (y_2_2 - 1) * 4 + (y_2_3 - 1) * 2 + y_2_4 ;
                                        
                                        hold_var = (y_1_1 - 1 ) * 32 + (y_1_2 - 1) * 16 + (y_1_3 - 1) * 8 + (y_2_1 - 1) * 4 + (y_2_2 - 1) * 2 + y_2_3;
                                        f_u_1_u_2_given_y_11_y_12_y_13_y_21_y_22_y_23 = f(: , : , hold_var) ;
                                        
                                        
                                        u_1_index = find (T_u_1(: , 5 + y_2_prime) == x_1_4 ) ;
                                        u_2_index = find (T_u_2(: , 5 + y_1_prime) == x_2_4) ;
                                        
                                        for u_1_i = 1 : length(u_1_index)
                                            u_1 = T_u_1(u_1_index(u_1_i) , 1) ;
                                            x_1_3 = T_u_1(u_1_index(u_1_i) , 1 + y_2_12) ;
                                            
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
                                                    * delta ...
                                                    * f_u_1_u_2_given_y_11_y_12_y_13_y_21_y_22_y_23(u_1_index(u_1_i) , u_2_index(u_2_i))' .* ((hold_var).^2 ...
                                                    + (u_2 - codebook_user_2(y_1 , : , u_1_index(u_1_i))).^2);
                                                
                                                summation_1 = summation_1 + P_bit_3(y_1_prime , y_2_prime) ...
                                                    * Pr_z_1i_given_z_1i_1(binary_z_13 + 1 , binary_z_14 + 1) ...
                                                    * Pr_z_2i_given_z_2i_1(binary_z_23 + 1 , binary_z_24 + 1) ...
                                                    * delta ...
                                                    * f_u_1_u_2_given_y_11_y_12_y_13_y_21_y_22_y_23(u_1_index(u_1_i) , u_2_index(u_2_i))' .* ((hold_var).^ 2) ;
                                                
                                                summation_2 = summation_2 + P_bit_3(y_1_prime , y_2_prime) ...
                                                    * Pr_z_1i_given_z_1i_1(binary_z_13 + 1 , binary_z_14 + 1) ...
                                                    * Pr_z_2i_given_z_2i_1(binary_z_23 + 1 , binary_z_24 + 1) ...
                                                    * delta ...
                                                    * f_u_1_u_2_given_y_11_y_12_y_13_y_21_y_22_y_23(u_1_index (u_1_i) , u_2_index(u_2_i))' .* ((u_2 - codebook_user_2(y_1 , : , u_1_index(u_1_i))).^2) ;
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end
D = summation ;
D_1 = summation_1 ;
D_2 = summation_2 ;
end