function P_bit_3 = P_y_11_y_12_y_13_y_21_y_22_y_23(P_bit_2 , Pr_z_1i_given_z_1i_1 , Pr_z_2i_given_z_2i_1 ,  T , f , delta)
summation = 0 ;
P_bit_3 = zeros (8 , 8) ;
for y_1_1 = 1 : 2
    for y_1_2 = 1 : 2
        y_1_prime = (y_1_1 - 1) * 2 + y_1_2 ;
        for y_1_3 = 1 : 2
            y_1 = (y_1_1 - 1) * 4 + (y_1_2 - 1) * 2 + y_1_3 ;
            for y_2_1 = 1 : 2
                for y_2_2 = 1 : 2
                    y_2_prime = (y_2_1 - 1 ) * 2 + y_2_2 ;
                    hold_var = (y_1_1 - 1) * 8 + (y_1_2 - 1) * 4 + (y_2_1 - 1) * 2  + y_2_2 ;
                    f_u_1_u_2_given_y_11_y_12_y_21_y_22 = f(: , : , hold_var) ;
                    for y_2_3 = 1 : 2
                        y_2 = (y_2_1 - 1) * 4 + (y_2_2 - 1) * 2 + y_2_3 ;
                        
                        for x_1_3 = 1 : 2
                            for x_2_3 = 1 : 2
                                
                                u_1_index = find (T (: , 3 + y_2_prime) == x_1_3) ;
                                u_2_index = find (T (: , 10 + y_1_prime) == x_2_3) ;
                                for u_1_i = 1 : length(u_1_index)
                                    x_1_2 = T(u_1_index(u_1_i) , 1 + y_2_1) ;
                                    for u_2_i = 1 : length(u_2_index)
                                        x_2_2 = T(u_2_index(u_2_i) , 8 + y_1_1) ;
                                        
                                        binary_z_13 = xor(x_2_3 - 1 , y_1_3 - 1 ) ;
                                        binary_z_12 = xor(x_2_2 - 1 , y_1_2 - 1) ;
                                        
                                        binary_z_23 = xor(x_1_3 - 1 , y_2_3 - 1) ;
                                        binary_z_22 = xor(x_1_2 - 1 , y_2_2 - 1) ;
                                        
                                        
                                        summation = summation + ...
                                            P_bit_2(y_1_prime , y_2_prime) ...
                                            * Pr_z_1i_given_z_1i_1(binary_z_12 + 1 , binary_z_13 + 1) ...
                                            * Pr_z_2i_given_z_2i_1(binary_z_22 + 1 , binary_z_23 + 1) ...
                                            * delta ...
                                            * f_u_1_u_2_given_y_11_y_12_y_21_y_22(u_1_index(u_1_i) , u_2_index(u_2_i))';
                                    end
                                end
                            end
                        end
                        P_bit_3(y_1 , y_2) = summation ;
                        summation = 0 ;
                    end
                end
            end
        end
    end
end
P_bit_3 = P_bit_3 ./ sum(sum(P_bit_3)) ;
end