function f_u_1_u_2_given_y_1_1_y_2_1 = generate_pdf_rate_2(y_1_1 , y_2_1 , T , Pr , delta , f)
%% numerator
numerator = zeros (length(T) , length(T)) ;
for u_1_index = 1 : length(T)
    for u_2_index = 1 : length(T)

        x_1_1 = T(u_1_index , 2) ;
        x_2_1 = T(u_2_index , 4) ;
        
        y = (y_1_1 - 1) * 2 + y_2_1 ;
        x = (x_2_1 - 1) * 2 + x_1_1 ;
        
        numerator(u_1_index , u_2_index) = Pr(x , y) * f(u_1_index , u_2_index) ;
        
    end
end
%% denominator
summation = 0 ;
for x_1_1 = 1 : 2
    for x_2_1 = 1 : 2
        y = (y_1_1 - 1) * 2 + y_2_1 ;
        x = (x_2_1 - 1) * 2 + x_1_1 ;
        
        u_1_index_x_1_1 = find (T(: , 2) == x_1_1) ;
        
        u_2_index_x_2_1 = find (T(: , 4) == x_2_1) ;
        for u_1_i = 1 : length(u_1_index_x_1_1)
        
            summation = summation + Pr (x , y) * sum(f(u_1_index_x_1_1(u_1_i) , u_2_index_x_2_1)') * delta ;
        end
    end
end
denominator = summation ;
f_u_1_u_2_given_y_1_1_y_2_1 = numerator ./ denominator ;
f_u_1_u_2_given_y_1_1_y_2_1 = f_u_1_u_2_given_y_1_1_y_2_1 ./ (sum(sum(f_u_1_u_2_given_y_1_1_y_2_1)) * delta) ;
end