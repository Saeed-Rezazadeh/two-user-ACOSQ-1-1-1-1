function f_u_1_u_2_given_y_11_y_12_y_13_y_21_y_22_y_23 = generate_pdf_rate_4(y_1_1 , y_1_2 , y_1_3 , y_2_1 , y_2_2 , y_2_3 , T , f , Pr_z_1i_given_z_1i_1 , Pr_z_2i_given_z_2i_1 , delta)
numerator = zeros (length(T) , length(T)) ;
y_1_prime = (y_1_1 - 1) * 2 + y_1_2 ;
y_2_prime = (y_2_1 - 1) * 2 + y_2_2 ;

for u_1_index = 1 : length(T)
    for u_2_index = 1 : length(T)
        
        
        
        
        x_1_3 = T(u_1_index , 3 + y_2_prime) ;
        x_1_2 = T(u_1_index , 1 + y_2_1) ;
        
        x_2_3 = T(u_2_index , 10 + y_1_prime) ;
        x_2_2 = T(u_2_index , 8 + y_1_1) ;
        
        binary_z_13 = xor(x_2_3 - 1 , y_1_3 - 1 ) ;
        binary_z_12 = xor(x_2_2 - 1 , y_1_2 - 1) ;
        
        binary_z_23 = xor(x_1_3 - 1 , y_2_3 - 1) ;
        binary_z_22 = xor(x_1_2 - 1 , y_2_2 - 1) ;
        
        numerator(u_1_index , u_2_index) = f(u_1_index , u_2_index) ...
            * Pr_z_1i_given_z_1i_1(binary_z_12 + 1 , binary_z_13 + 1) ...
            * Pr_z_2i_given_z_2i_1(binary_z_22 + 1 , binary_z_23 + 1);
    end
end

summation = 0 ;
for x_1_3 = 1 : 2
    for x_2_3 = 1 : 2
        
        
        binary_z_13 = xor(x_2_3 - 1 , y_1_3 - 1 ) ;
        binary_z_12 = xor(x_2_2 - 1 , y_1_2 - 1) ;
        
        binary_z_23 = xor(x_1_3 - 1 , y_2_3 - 1) ;
        binary_z_22 = xor(x_1_2 - 1 , y_2_2 - 1) ;
        
        
        u_1_index = find (T(: , 3 + y_2_prime) == x_1_3) ;
        u_2_index = find (T(: , 10 + y_1_prime) == x_2_3) ;
        for u_1_i = 1 : length(u_1_index)
            for u_2_i = 1 : length(u_2_index)
                summation = summation ...
                    + Pr_z_1i_given_z_1i_1(binary_z_12 + 1 , binary_z_13 + 1)...
                    * Pr_z_2i_given_z_2i_1(binary_z_22 + 1 , binary_z_23 + 1)...
                    * delta * f(u_1_index(u_1_i) , u_2_index(u_2_i)) ;
            end
        end
    end
end
denominator = summation ;
f_u_1_u_2_given_y_11_y_12_y_13_y_21_y_22_y_23 = numerator ./ denominator ;
f_u_1_u_2_given_y_11_y_12_y_13_y_21_y_22_y_23 = f_u_1_u_2_given_y_11_y_12_y_13_y_21_y_22_y_23 ./ (sum(sum(f_u_1_u_2_given_y_11_y_12_y_13_y_21_y_22_y_23)) * delta) ;
end