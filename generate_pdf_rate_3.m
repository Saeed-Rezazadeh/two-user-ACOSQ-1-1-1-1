function f_u_1_u_2_given_y_11_y_12_y_21_y_22 = generate_pdf_rate_3(y_1_1 , y_1_2 , y_2_1 , y_2_2 , T , f , Pr_z_1i_given_z_1i_1 , Pr_z_2i_given_z_2i_1  , delta)
%% numerator
numerator = zeros (length(T) , length(T)) ;
for u_1_index = 1 : length(T)
    for u_2_index = 1 : length(T)
        x_1_2 = T(u_1_index , 2 + y_2_1) ;
        x_2_2 = T(u_2_index , 6 + y_1_1) ;
        
        x_1_1 = T(u_1_index , 2 ) ;
        x_2_1 = T(u_2_index , 6 ) ;
        
        binary_z_11 = xor(x_2_1 - 1 , y_1_1 - 1 ) ;
        binary_z_12 = xor(x_2_2 - 1 , y_1_2 - 1) ;
        
        binary_z_21 = xor(x_1_1 - 1 , y_2_1 - 1) ;
        binary_z_22 = xor(x_1_2 - 1 , y_2_2 - 1) ;
        
        
        numerator(u_1_index , u_2_index) = f(u_1_index , u_2_index) ...
            * Pr_z_1i_given_z_1i_1(binary_z_11 + 1 , binary_z_12 + 1) ...
            * Pr_z_2i_given_z_2i_1(binary_z_21 + 1 , binary_z_22 + 1);
    end
end

%% denominator
summation = 0 ;
for x_1_2 = 1 : 2
    for x_2_2 = 1 : 2

        u_1_index_x_1_2 = find (T(: , 2 + y_2_1) == x_1_2) ;
        u_2_index_x_2_2 = find (T(: , 6 + y_1_1) == x_2_2) ;
        for u_1_i = 1 : length(u_1_index_x_1_2)
            for u_2_i = 1 : length(u_2_index_x_2_2)
                
                x_1_1 = T(u_1_index , 2 ) ;
                x_2_1 = T(u_2_index , 6 ) ;
                
                binary_z_11 = xor(x_2_1 - 1 , y_1_1 - 1 ) ;
                binary_z_12 = xor(x_2_2 - 1 , y_1_2 - 1) ;
                
                binary_z_21 = xor(x_1_1 - 1 , y_2_1 - 1) ;
                binary_z_22 = xor(x_1_2 - 1 , y_2_2 - 1) ;
                
                summation = summation ...
                    + Pr_z_1i_given_z_1i_1(binary_z_11 + 1 , binary_z_12 + 1) ...
                    * Pr_z_2i_given_z_2i_1(binary_z_21 + 1 , binary_z_22 + 1)...
                    * delta * f(u_1_index_x_1_2(u_1_i) , u_2_index_x_2_2(u_2_i)) ;
            end
        end
    end
end
denominator = summation ;
f_u_1_u_2_given_y_11_y_12_y_21_y_22 = numerator ./ denominator ;
f_u_1_u_2_given_y_11_y_12_y_21_y_22 = f_u_1_u_2_given_y_11_y_12_y_21_y_22 ./ (sum(sum(f_u_1_u_2_given_y_11_y_12_y_21_y_22)) * delta ) ;
end