function [D , D_1 , D_2] = distortion_2(f , codebook_user_1 , codebook_user_2 , Pr_z_1i_given_z_1i_1 , Pr_z_2i_given_z_2i_1 , P_bit_1 , T_u_1 ,T_u_2 , delta)
summation = 0 ;
summation_1 = 0 ;
summation_2 = 0 ;
parfor y_1_1 = 1 : 2
    for y_2_1 = 1 : 2
         hold_var = (y_1_1 - 1 ) * 2 + y_2_1 ;
         f_u_1_u_2_given_y_1_1_y_2_1 = f (: , : , hold_var) ;
        for x_1_2 = 1 : 2
            for x_2_2 = 1 : 2
                for y_1_2 = 1 : 2
                    y_1 = (y_1_1 - 1) * 2 + y_1_2 ;
                    for y_2_2 = 1 : 2
                        y_2 = (y_2_1 - 1) * 2 + y_2_2 ;
                        u_1_index = find (T_u_1(: , 2 + y_2_1) == x_1_2) ;
                        u_2_index = find (T_u_2(: , 2 + y_1_1) == x_2_2) ;
                        for u_1_i = 1 : length(u_1_index)
                            x_1_1 = T_u_1(u_1_index(u_1_i) , 2) ;
                            u_1 = T_u_1(u_1_index(u_1_i) , 1) ;
                            for u_2_i = 1 : length(u_2_index)
                                u_2 = T_u_2(u_2_index(u_2_i) , 1) ;
                                x_2_1 = T_u_2(u_2_index(u_2_i) , 2) ;
                                binary_z_11 = xor(x_2_1 - 1 , y_1_1 - 1 ) ;
                                binary_z_12 = xor(x_2_2 - 1 , y_1_2 - 1) ;
                                
                                binary_z_21 = xor(x_1_1 - 1 , y_2_1 - 1) ;
                                binary_z_22 = xor(x_1_2 - 1 , y_2_2 - 1) ;
                                
                                
                                
                                hold_var_1 = u_1 - codebook_user_1(y_2 , u_2_index(u_2_i)) ;
                                hold_var_1 = hold_var_1(:) ;
                                
                                summation = summation  + P_bit_1(y_1_1 , y_2_1) ...
                                    * Pr_z_1i_given_z_1i_1(binary_z_11 + 1 , binary_z_12 + 1) ...
                                    * Pr_z_2i_given_z_2i_1(binary_z_21 + 1 , binary_z_22 + 1) ...
                                    * delta * f_u_1_u_2_given_y_1_1_y_2_1(u_1_index(u_1_i) , u_2_index(u_2_i))' ...
                                    .* ((hold_var_1).^2 + (u_2 - codebook_user_2(y_1 , : , u_1_index(u_1_i))).^2) ;
                                
                                summation_1 = summation_1 + P_bit_1(y_1_1 , y_2_1) ...
                                    * Pr_z_1i_given_z_1i_1(binary_z_11 + 1 , binary_z_12 + 1) ...
                                    * Pr_z_2i_given_z_2i_1(binary_z_21 + 1 , binary_z_22 + 1) ...
                                    * delta * f_u_1_u_2_given_y_1_1_y_2_1(u_1_index(u_1_i) , u_2_index(u_2_i))' ...
                                    .* (hold_var_1).^2 ;
                                
                                summation_2 = summation_2 + P_bit_1(y_1_1 , y_2_1) ...
                                    * Pr_z_1i_given_z_1i_1(binary_z_11 + 1 , binary_z_12 + 1) ...
                                    * Pr_z_2i_given_z_2i_1(binary_z_21 + 1 , binary_z_22 + 1) ...
                                    * delta * f_u_1_u_2_given_y_1_1_y_2_1(u_1_index(u_1_i) , u_2_index(u_2_i))' ...
                                    .*(u_2 - codebook_user_2(y_1 , : , u_1_index(u_1_i))).^2 ;
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