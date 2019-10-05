function p_bit_2 = P_y_11_y_12_y_21_y_22 (P_bit_1 , Pr_z_1i_given_z_1i_1 , Pr_z_2i_given_z_2i_1 , T , f , delta)
summation = 0 ;
p_bit_2 = zeros (4 , 4) ;
for y_1_1 = 1 : 2
    for y_1_2 = 1 : 2
        y_1 = (y_1_1 - 1) * 2  + y_1_2 ;
        for y_2_1 = 1 : 2
            y = (y_1_1 - 1) * 2 + y_2_1 ;
            f_u_1_u_2_given_y_11_y_21 = f(: , : , y) ;
            for y_2_2 = 1 : 2
                y_2 = (y_2_1 - 1) * 2 + y_2_2 ;
                for x_1_2 = 1 : 2
                    for x_2_2 = 1 : 2
                        u_1_index = find (T(: , 2 + y_2_1) == x_1_2) ;
                        u_2_index = find (T(: , 6 + y_1_1) == x_2_2) ;
                        
                        
                        
                        for u_1_i = 1 : length(u_1_index)
                            for u_2_i = 1 : length(u_2_index)
                                
                                x_1_1 = T(u_1_index(u_1_i) , 2 ) ;
                                x_2_1 = T(u_2_index(u_2_i) , 6 ) ;
                                binary_z_11 = xor(x_2_1 - 1 , y_1_1 - 1 ) ;
                                binary_z_12 = xor(x_2_2 - 1 , y_1_2 - 1) ;
                                
                                binary_z_21 = xor(x_1_1 - 1 , y_2_1 - 1) ;
                                binary_z_22 = xor(x_1_2 - 1 , y_2_2 - 1) ;
                                
                                
                                summation = summation  + P_bit_1(y_1_1 , y_2_1) ...
                                    * Pr_z_1i_given_z_1i_1(binary_z_11 + 1 , binary_z_12 + 1) ...
                                    * Pr_z_2i_given_z_2i_1(binary_z_21 + 1 , binary_z_22 + 1) ...
                                    * delta * f_u_1_u_2_given_y_11_y_21(u_1_index(u_1_i) , u_2_index(u_2_i)) ;
                            end
                        end
                    end
                end
                p_bit_2(y_1 , y_2) = summation ;
                summation = 0 ;
            end
        end
    end
end
p_bit_2 = p_bit_2 ./ sum(sum(p_bit_2));
end