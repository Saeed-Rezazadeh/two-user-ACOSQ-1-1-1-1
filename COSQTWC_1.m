function [SDR , SDR_1 , SDR_2 , T , codebook_user_1 , codebook_user_2] = COSQTWC_1 (p , Pr , T , width , f)
codebook_user_1 = zeros (2 , 1 , length(T )) ;
codebook_user_2 = zeros (2 , 1 , length(T )) ;
delta = width(2) * width(1) ;

%% pdfs
f_u_1_given_u_2 = zeros(length(T) , length(T)) ;
u_1 = T(: , 1) ;
for u_2_index = 1 : length(T)
    u_2 = T(u_2_index , 3) ;
    f_u_1_given_u_2(: , u_2_index) = 1 ./ (sqrt (2 .* pi .* (1 - p .^ 2))) .* exp (- 1 ./ (2 .* (1 - p .^ 2)) .* (u_1 - p .* u_2) .^ 2) ;
    f_u_1_given_u_2(: , u_2_index) = f_u_1_given_u_2(: , u_2_index) ./ (sum(f_u_1_given_u_2(: , u_2_index)) * width(1)) ;
end

f_u_2_given_u_1 = zeros(length(T) , length(T)) ;
u_2 = T(: , 3) ;
for u_1_index = 1 : length(T)
    u_1 = T(u_1_index , 1) ;
    f_u_2_given_u_1(u_1_index , :) = 1 ./ (sqrt (2 .* pi .* (1 - p .^ 2))) .* exp (- 1 ./ (2 .* (1 - p .^ 2)) .* (u_2 - p .* u_1) .^ 2) ;
    f_u_2_given_u_1(u_1_index , :) = f_u_2_given_u_1(u_1_index , :) ./ (sum(f_u_2_given_u_1(u_1_index , :)) * width(2)) ;
end

Threshold = 0.001 ;
D = [1 2] ;
FileID = fopen ('Results.txt', 'a') ;

while abs(D(2) - D(1)) / D(2) > (Threshold)
    D (1) = D(2) ;
    %% Optimal Codebook for user 1
    parfor u_2_index =1 : length(T)
        u_2 = T(u_2_index , 3) ;
        i_2 = T(u_2_index , 4) ;
        for j_2 = 1 : 2
            numerator = 0 ;
            denominator = 0 ;
            for i_1 = 1 : 2
                for j_1 = 1 : 2
                    i = (i_2 - 1) * 2 + i_1 ;
                    j = (j_1 - 1) * 2 + j_2 ;
                    
                    u_1_index = find (T(: , 2) == i_1) ;
                    u_1 = T(u_1_index , 1) ;
                    
                    
                    numerator = numerator + Pr(i , j) * sum(u_1 .* f_u_1_given_u_2(u_1_index , u_2_index)) ;
                    denominator = denominator + Pr(i , j) * sum(f_u_1_given_u_2(u_1_index , u_2_index)) ;
                end
            end
            if (numerator == 0 && denominator ==0 )
                codebook_user_1(j_2 , : , u_2_index) = p * u_2 ;
            else
                codebook_user_1(j_2 , : , u_2_index) = numerator / denominator ;
            end
        end
    end
    %% Optimal Codebook for user 2
    parfor u_1_index = 1 : length(T)
        u_1 = T(u_1_index , 1) ;
        i_1 = T(u_1_index , 2) ;
        for j_1 = 1 : 2
            numerator = 0 ;
            denominator = 0 ;
            for i_2 = 1 : 2
                for j_2 = 1 : 2
                    i = (i_2 - 1) * 2 + i_1 ;
                    j = (j_1 - 1) * 2 + j_2 ;
                    
                    u_2_index = find (T(: , 4) == i_2) ;
                    u_2 = T(u_2_index , 3) ;
                    
                    
                    numerator = numerator + Pr(i , j) * sum(u_2 .* f_u_2_given_u_1(u_1_index , u_2_index)') ;
                    denominator = denominator + Pr(i , j) * sum(f_u_2_given_u_1(u_1_index , u_2_index )') ;
                end
            end
            if(numerator == 0 && denominator == 0)
                codebook_user_2(j_1 , : , u_1_index) = p * u_1 ;
            else
                codebook_user_2(j_1 , : , u_1_index) = numerator / denominator ;
            end
        end
    end
    
    %% Optimal Partitions for user 1
    parfor u_1_index = 1 : length(T)
        u_1 = T(u_1_index , 1) ;
        summation = 0 ;
        temp = zeros (2 , 1) ;
        for i_1 = 1 : 2
            for i_2 = 1 : 2
                for j_1 = 1 : 2
                    for j_2 = 1 : 2
                        j = (j_1 - 1) * 2 + j_2 ;
                        i = (i_2 - 1) * 2 + i_1 ;
                        
                        u_2_index = find (T(: , 4) == i_2) ;
                        u_2 = T(u_2_index , 3) ;
                        
                        hold_var = u_1 - codebook_user_1(j_2 , : , u_2_index) ;
                        hold_var = hold_var(:) ;
                        
                        
                        summation = summation + Pr(i , j) * width(2) * sum(f_u_2_given_u_1(u_1_index , u_2_index)' .* ...
                            ((hold_var) .^ 2 + (u_2 - codebook_user_2(j_1 , : , u_1_index)) .^ 2)) ;
                        
                    end
                end
            end
            temp (i_1) = summation ;
            summation = 0 ;
        end
        [~ , partition_index] = min(temp) ;
        T_u_1(u_1_index , 1) = partition_index ;
    end
    T(: , 2) = T_u_1 (: , 1);
    
    %% Optimal Partition user 2
    parfor u_2_index = 1 : length(T)
        summation = 0 ;
        temp = zeros(2 , 1) ;
        u_2 = T(u_2_index , 3) ;
        for i_2 = 1 : 2
            for i_1 = 1 : 2
                for j_1 = 1 : 2
                    for j_2 = 1 : 2
                        i = (i_2 - 1) * 2 + i_1 ;
                        j = (j_1 - 1) * 2 + j_2 ;
                        
                        u_1_index = find (T(: , 2) == i_1) ;
                        u_1 = T(u_1_index , 1) ;
                        
                        hold_var = u_2 - codebook_user_2(j_1 , : , u_1_index) ;
                        hold_var = hold_var(:) ;
                        
                        
                        summation = summation + Pr(i , j) * width(1) * sum(f_u_1_given_u_2(u_1_index , u_2_index) ...
                            .* ((u_1 - codebook_user_1(j_2 , : , u_2_index)) .^ 2 + (hold_var) .^ 2)) ;
                        
                    end
                end
            end
            temp (i_2) = summation ;
            summation = 0 ;
        end
        [~ , partition_index ] = min(temp) ;
        T_u_2(u_2_index , 1) = partition_index ;
    end
    T(: , 4) = T_u_2(: , 1) ;
    
    %% Distortion
    [D(2) , D_1 , D_2] = distortion_1 (Pr , codebook_user_1 , codebook_user_2 , T , delta , f) ;
    fprintf (FileID , '\n Overall: \n') ;
    fprintf (FileID , ' %f ' , D(2)) ;
    fprintf (FileID , ' %f ' , D_1) ;
    fprintf (FileID , ' %f ' , D_2) ;
end
Distortion = D(2) ;
SDR = 10 * log10(2 / Distortion) ;
SDR_1 = 10 * log10(1 / D_1) ;
SDR_2 = 10 * log10(1 / D_2) ;
fclose (FileID) ;
end