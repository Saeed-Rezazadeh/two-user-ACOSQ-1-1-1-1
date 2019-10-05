%% This is the MATLAB script ACOSQ with r_1 = r_2 = (1 1 1 1)design extended to a two-user system. The simulations results are provided in Appendix A.
clc;
close all
clear all
%% The Results.txt
% This file contains the ultimate SDR values for different channel parameters delta and epsilon.
% Also, since the proposed ACOSQ is an iterative algorithm, the distirtion
% value for every iteration is provided in this file for given a epsilon and delta.
FileID = fopen ('Results.txt', 'a') ;



%% Sources Correlation
rho = [0.5 , 0.9] ;

%% Noise parameters
epsilon_1 = unique ([ 10^-3 , 0.005 0.01  0.05  0.1]);
epsilon_2 = unique ([ 10^-3 , 0.005 0.01  0.05  0.1]);
SIZE = length(epsilon_1) ;
i_index = [1 : SIZE , SIZE : -1 : 1  , 1 : SIZE , SIZE : -1 : 1] ;

Delta = [0 10] ;
% Channel MODE
MODE = 1 ; % 1 = Orthogonal  2 = Additive 3 = Multiplicative


%% Initialize parameters
final_SDR_rate_4 = zeros (length(Delta) , length(rho) , length(i_index)) ;
final_SDR_rate_3 = zeros (length(Delta) , length(rho) , length(i_index)) ;
final_SDR_rate_2 = zeros (length(Delta) , length(rho) , length(i_index)) ;
final_SDR_rate_1 = zeros (length(Delta) , length(rho) , length(i_index)) ;

SDR_rate_1 = zeros (length(Delta) , length(rho) , length(i_index)) ;
SDR_1_rate_1 = zeros (length(Delta) , length(rho) , length(i_index)) ;
SDR_2_rate_1 = zeros (length(Delta) , length(rho) , length(i_index)) ;

SDR_rate_2 = zeros (length(Delta) , length(rho) , length(i_index)) ;
SDR_1_rate_2 = zeros (length(Delta) , length(rho) , length(i_index)) ;
SDR_2_rate_2 = zeros (length(Delta) , length(rho) , length(i_index)) ;

SDR_rate_3 = zeros (length(Delta) , length(rho) , length(i_index)) ;
SDR_1_rate_3 = zeros (length(Delta) , length(rho) , length(i_index)) ;
SDR_2_rate_3 = zeros (length(Delta) , length(rho) , length(i_index)) ;

SDR_rate_4 = zeros (length(Delta) , length(rho) , length(i_index)) ;
SDR_1_rate_4 = zeros (length(Delta) , length(rho) , length(i_index)) ;
SDR_2_rate_4 = zeros (length(Delta) , length(rho) , length(i_index)) ;

% Determine the correlation between the two sources.
for j = 1 : length(rho)
    p = rho (j) ;
    % Determine the amount of noise correlation.
    for h = 1 : length(Delta)
        delta = Delta(h) ;
        % Determine the noise parameters in the either of direction of
        % transmissions.
        for k = 1 : length(i_index)
            i = i_index(k) ;
            
            %% Load Data
            % Initial partitions are obtained from half-duplex ACOSQ for one-way
            % system with side information at the decoder. Then the
            % half-duplex design is used to initialize the two-user ACOSQ
            numLevel = 2 ;
            Data = ['Inputs\T_k_' num2str(k) '_delta_' num2str(delta) '_p_0_' num2str(10 * p)] ;
            load (Data) ;
            
            width(1) = T(2 , 1) - T (1 , 1) ;
            width(2) = width(1) ;
            delta_u = width(1) * width(2) ;
            joint_T = cat(2 , T(: , [1 2]) , T(: , [1 2])) ;
            f = zeros(length(joint_T) , length(joint_T)) ;
            u_1 = joint_T(: , 1) ;
            
            % Compute the source pdf f_u_1_u_2.
            for u_2_index = 1 : length(joint_T)
                u_2 = joint_T(u_2_index , 3) ;
                f(: , u_2_index ) = 1 ./ (2 .* pi .* sqrt (1 - p .^ 2)) .* exp (-1 ./ (2 .* (1 - p .^ 2)) .* (u_1 .^ 2 + u_2 .^ 2 - 2 .* p .* u_1 .* u_2)) ;
                
            end
            f = f ./ (sum(sum(f)) * delta_u) ;
            
            %% Generate Channel transission distributions based on the Polya Contagion Urn process.
            Pr_z_1i_given_z_1i_1 = [(1 - epsilon_1(i) + delta) / (1 + delta)  , epsilon_1(i) / (1 + delta) ;
                (1 - epsilon_1(i)) / (1 + delta) , (epsilon_1(i) + delta ) / (1 + delta)] ;
            
            Pr_z_2i_given_z_2i_1 = [(1 - epsilon_2(i) + delta) / (1 + delta)  , epsilon_2(i) / (1 + delta) ;
                (1 - epsilon_2(i)) / (1 + delta) , (epsilon_2(i) + delta ) / (1 + delta)] ;
            
            
            [Probability] = channel_statistics (MODE  , epsilon_1 , epsilon_2 , numLevel , i , delta) ;
            
            %% COSQTWC for step #1
            % Design a 1-bit COSQ for the sources with joint pdf computed in line
            % 78.
            [SDR_rate_1(h , j , k) , SDR_1_rate_1(h , j , k) , SDR_2_rate_1(h , j , k) , joint_T_rate_1 , codebook_user_1 , codebook_user_2] =...
                COSQTWC_1 (p , Probability , joint_T , width , f) ;
            fprintf (FileID , '\n i = %4.2f \n' , i ) ;
            fprintf (FileID , '\n SDR values: \n') ;
            fprintf (FileID , ' SDR = %7.4f SDR1 = %7.4f SDR2 = %7.4f' , SDR_rate_1(h , j , k) , SDR_1_rate_1(h , j , k) , SDR_2_rate_1(h , j , k)) ;
            fprintf (FileID , '\n rate 1 \n') ;
            
            
            
            %% COSQTWC for step #2
            % Compute the probabilities Pr(Y_11Y_21 = y_11y_21) for
            % y_11y_21 = 00, 01 , 10 , and 11.
            P_bit_1 = P_y_11_y_21(Probability , joint_T_rate_1 , width , f) ;
            
            % Compute the conditional source pdf f_u_1_u_2_given_y_11_y_21
            % by extending (4.8) to the two-user system setup described in
            % Appendix A.
            f_u_1_u_2_given_y_11_y_21 = zeros (length(joint_T_rate_1) , length(joint_T_rate_1) , 4);
            for y_1_1 = 1 : 2
                for y_2_1 = 1 : 2
                    y_1 = (y_1_1 - 1 ) * 2 + y_2_1 ;
                    f_u_1_u_2_given_y_11_y_21(: , : , y_1) = ...
                        generate_pdf_rate_2(y_1_1 , y_2_1 , joint_T_rate_1 , Probability , delta_u  , f) ;
                end
            end
            
            numLevel = 4 ;
            T_u_1 = T(: , [1 2 3 4]) ;
            T_u_2 = T(: , [1 2 3 4]) ;
            % Design a 1-bit COSQ for the joint source pdf computed in Line
            % 117.
            [SDR_rate_2(h , j , k)  , SDR_1_rate_2(h , j , k)  , SDR_2_rate_2(h , j , k)  , joint_T_rate_2 , codebook_user_1 , codebook_user_2] ...
                = COSQTWC_2 (f_u_1_u_2_given_y_11_y_21 , p , ...
                T_u_1 , T_u_2 , Pr_z_1i_given_z_1i_1 , Pr_z_2i_given_z_2i_1 , ...
                P_bit_1 , width ,  codebook_user_1 , codebook_user_2) ;
            fprintf (FileID , '\n i = %4.2f \n' , i ) ;
            fprintf (FileID , '\n SDR values: \n') ;
            fprintf (FileID , ' SDR = %7.4f SDR1 = %7.4f SDR2 = %7.4f' , SDR_rate_2(h , j , k) , SDR_1_rate_2(h , j , k) , SDR_2_rate_2(h , j , k)) ;
            fprintf (FileID , '\n rate 2 \n') ;
            
            %% COSQ for rate 3
            % Compute the probabilities Pr(Y_11Y_12Y_21Y_22 = y_11y_12y_21y_22) for
            % y_11y_12y_21y_22 = 0000 , 0001 , ... , and 1111.
            P_bit_2 = P_y_11_y_12_y_21_y_22 (P_bit_1 , Pr_z_1i_given_z_1i_1 , Pr_z_2i_given_z_2i_1 ...
                ,  joint_T_rate_2 , f_u_1_u_2_given_y_11_y_21 , width(2) * width(1)) ;
            
            % Compute the conditional source pdf f_u_1_u_2_given_y_11_y_12_y_21_y_22
            % by extending (4.8) to the two-user system setup described in
            % Appendix A.
            f_u_1_u_2_given_y_11_y_12_y_21_y_22 = zeros (length(joint_T_rate_2) , length(joint_T_rate_2) , 16);
            for y_1_1 = 1 : 2
                for y_1_2 = 1 : 2
                    for y_2_1 = 1 : 2
                        y_1 = (y_1_1 - 1) * 2 + y_2_1 ;
                        for y_2_2 = 1 : 2
                            y_2 = (y_1_1 - 1) * 8 + (y_1_2 - 1) * 4 + (y_2_1 - 1) * 2 + y_2_2 ;
                            pdf_given_y_11_y_21 = f_u_1_u_2_given_y_11_y_21(: , : , y_1) ;
                            hold_var = ...
                                generate_pdf_rate_3(y_1_1 , y_1_2 , y_2_1 , y_2_2 ...
                                , joint_T_rate_2 , pdf_given_y_11_y_21 ...
                                , Pr_z_1i_given_z_1i_1 , Pr_z_2i_given_z_2i_1 , width(2) * width(1));
                            
                            f_u_1_u_2_given_y_11_y_12_y_21_y_22(: , : , y_2) = hold_var ;
                        end
                    end
                end
            end
            T_u_1 = T(: , [1 3 4 5 6 7 8]) ;
            T_u_2 = T(: , [1 3 4 5 6 7 8]) ;
            % Design a 1-bit COSQ for the joint source pdf computed in Line
            % 158.
            [SDR_rate_3(h , j , k) , SDR_1_rate_3(h , j , k) , SDR_2_rate_3(h , j , k) , joint_T_rate_3 , codebook_user_1 , codebook_user_2] = ...
                COSQTWC_3(p , f_u_1_u_2_given_y_11_y_12_y_21_y_22 , T_u_1 , T_u_2...
                , Pr_z_1i_given_z_1i_1 , Pr_z_2i_given_z_2i_1 , P_bit_2 ...
                , width ,  codebook_user_1 , codebook_user_2) ;
            fprintf (FileID , '\n i = %4.2f \n' , i ) ;
            fprintf (FileID , '\n SDR values: \n') ;
            fprintf (FileID , ' SDR = %7.4f SDR1 = %7.4f SDR2 = %7.4f' , SDR_rate_3(h , j , k) , SDR_1_rate_3(h , j , k) , SDR_2_rate_3(h , j , k)) ;
            fprintf (FileID , '\n rate 3 \n') ;
            
            
            %% COSQ for rate 4
            % Compute the probabilities Pr(Y_11Y_12Y_13Y_21Y_22Y_23 = y_11y_12y_13y_21y_22y_23) for
            % y_11y_12y_13y_21y_22y_23 = 000000 , 000001 , ... , and 111111.
            P_bit_3 = P_y_11_y_12_y_13_y_21_y_22_y_23 (P_bit_2 , Pr_z_1i_given_z_1i_1 , Pr_z_2i_given_z_2i_1 ,  joint_T_rate_3 , f_u_1_u_2_given_y_11_y_12_y_21_y_22 , width(2) * width(1)) ;
            
            % Compute the conditional source pdf f_u_1_u_2_y_11_y_12_y_13_y_21_y_22_y_23
            % by extending (4.8) to the two-user system setup described in
            % Appendix A.
            f_u_1_u_2_y_11_y_12_y_13_y_21_y_22_y_23 = zeros (length(joint_T_rate_3) , length(joint_T_rate_3) , 64) ;
            for y_1_1 = 1 : 2
                for y_1_2 = 1 : 2
                    for y_1_3 = 1 : 2
                        for y_2_1 = 1 : 2
                            for y_2_2 = 1 : 2
                                y_2 = (y_1_1 - 1) * 8 + (y_1_2 - 1) * 4 + (y_2_1 - 1) * 2 + y_2_2 ;
                                for y_2_3 = 1 : 2
                                    y_3 = (y_1_1 - 1 ) * 32 + (y_1_2 - 1) * 16 + (y_1_3 - 1) * 8 + (y_2_1 - 1) * 4 + (y_2_2 - 1) * 2 + y_2_3;
                                    pdf_given_y_11_y_12_y_21_y_22 = f_u_1_u_2_given_y_11_y_12_y_21_y_22(: , : , y_2) ;
                                    hold_var = ...
                                        generate_pdf_rate_4(y_1_1 , y_1_2 , y_1_3 , y_2_1 , y_2_2 , y_2_3 ...
                                        , joint_T_rate_3 , pdf_given_y_11_y_12_y_21_y_22 ...
                                        , Pr_z_1i_given_z_1i_1 , Pr_z_2i_given_z_2i_1 , width(2) * width(1));
                                    f_u_1_u_2_y_11_y_12_y_13_y_21_y_22_y_23(: , : , y_3) = hold_var;
                                end
                            end
                        end
                    end
                end
            end
            
            T_u_1 = T(: , [1 5 6 7 8 9 10 11 12 13 14 15 16]) ;
            T_u_2 = T(: , [1 5 6 7 8 9 10 11 12 13 14 15 16]) ;
            % Design a 1-bit COSQ for the joint source pdf computed in Line
            % 199.
            [SDR_rate_4(h , j , k) , SDR_1_rate_4(h , j , k) , SDR_2_rate_4(h , j , k) , joint_T_rate_4 , codebook_user_1 , codebook_user_2] = ...
                COSQTWC_4(p , f_u_1_u_2_y_11_y_12_y_13_y_21_y_22_y_23 ...
                , T_u_1 , T_u_2 , Pr_z_1i_given_z_1i_1 , Pr_z_2i_given_z_2i_1  ...
                , P_bit_3 , width ,  codebook_user_1 , codebook_user_2) ;
            fprintf (FileID , '\n i = %4.2f \n' , i ) ;
            fprintf (FileID , '\n SDR values: \n') ;
            fprintf (FileID , ' SDR = %7.4f SDR1 = %7.4f SDR2 = %7.4f' , SDR_rate_4(h , j , k) , SDR_1_rate_4(h , j , k) , SDR_2_rate_4(h , j , k)) ;
            fprintf (FileID , '\n rate 4 \n') ;
            
            % Save the ultimate quantizers for experimental results.
            Data = ['T\Data_k_' num2str(k) '_p_' num2str(j) '_noise_Delta_' num2str(delta)] ;
            save (Data , 'joint_T_rate_4' , 'joint_T_rate_3' , 'joint_T_rate_2' , 'joint_T_rate_1') ;
        end
        
        % Pick the best SDR values
        fprintf (FileID , '\nOver All SDR for rate 1\n') ;
        clear D_1 ;
        D_1 = zeros(SIZE , 1) ;
        for i = 1 : SIZE
            index = find (i_index == i) ;
            hold_var = SDR_rate_1(h , j , index) ;
            hold_var = hold_var(:) ;
            final_SDR_rate_1(h , j , i) = max(hold_var) ;
            fprintf (FileID , '\ni = %d' , i) ;
            D_1(i) = 10 ^(- final_SDR_rate_1(h , j , i) / 10) ;
            fprintf (FileID , '\nD = %f' , D_1(i)) ;
            fprintf (FileID , '\nSDR for rate 1 = %7.4f\n' , final_SDR_rate_1(h , j , i)) ;
        end
        fprintf (FileID , '\nOverall SDR for rate 2\n') ;
        clear D_2 ;
        D_2 = zeros(SIZE , 1) ;
        for i = 1 : SIZE
            index = find (i_index == i) ;
            hold_var = SDR_rate_2(h , j , index) ;
            hold_var = hold_var(:) ;
            final_SDR_rate_2(h , j , i) = max(hold_var) ;
            D_2(i) = 10 ^(- final_SDR_rate_2(h , j , i) / 10) ;
            fprintf (FileID , '\nD = %f' , D_2(i)) ;
            fprintf (FileID , '\nSDR for rate 2 = %7.4f\n' , final_SDR_rate_2(h , j , i)) ;
        end
        fprintf (FileID , '\nOverall SDR for rate 3\n') ;
        clear D_3 ;
        D_3 = zeros(SIZE , 1) ;
        for i = 1 : SIZE
            index = find (i_index == i) ;
            hold_var = SDR_rate_3(h , j , index) ;
            hold_var = hold_var(:) ;
            final_SDR_rate_3(h , j , i) = max(hold_var) ;
            D_3(i) = 10 ^(- final_SDR_rate_3(h , j , i) / 10) ;
            fprintf (FileID , '\nD = %f' , D_3(i)) ;
            fprintf (FileID , '\nSDR for rate 3 = %7.4f\n' , final_SDR_rate_3(h , j , i)) ;
        end
        
        fprintf (FileID , '\nOverall SDR for rate 4\n') ;
        clear D_4 ;
        D_4 = zeros(SIZE , 1) ;
        for i = 1 : SIZE
            index = find (i_index == i) ;
            hold_var = SDR_rate_4(h , j , index) ;
            hold_var = hold_var(:) ;
            final_SDR_rate_4(h , j , i) = max(hold_var) ;
            D_4(i) = 10 ^(- final_SDR_rate_4(h , j , i) / 10) ;
            fprintf (FileID , '\nD = %f' , D_4(i)) ;
            fprintf (FileID , '\nSDR for rate 4 = %7.4f\n' , final_SDR_rate_4(h , j , i)) ;
        end
        fprintf (FileID , '\n noise_Delta = %4.2f \n' , delta ) ;
    end
    fprintf (FileID , '\n p = %4.2f \n' , p ) ;
end

