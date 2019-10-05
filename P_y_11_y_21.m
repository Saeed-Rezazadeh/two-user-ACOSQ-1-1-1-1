function P = P_y_11_y_21(Pr , T , width , f)
delta = width(2) * width(1) ;
P = zeros (2 , 2) ;
summation = 0 ;
for y_1_1 = 1 : 2
    for y_2_1 = 1 : 2
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
        P(y_1_1 , y_2_1) = summation ;
        summation = 0 ;
    end
end
P = P ./ sum(sum(P));
end