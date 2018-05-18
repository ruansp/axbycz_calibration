% This is the main file for probabilistic AXB = YCZ problem: 
%  Comparisons between different algorithms with data scrambled
%
% Authors: Qianli Ma, qianli.ma622@gmail.com; 
%          Zachariah Goh, zach_goh@yahoo.com
% Modifications: Sipu Ruan, ruansp@jhu.edu
clear; close all; clc

%% Add file dependencies
addpath util/
addpath solvers/

%%
counter = 0;
gmean = [0; 0; 0; 0; 0 ;0];
Cov = eye(6,6);
k = 0.1;
num = 1; % number of simulations
Num = 80; % number of data
rate = 0:20:100;

%% generate random X, Y and Z
opt_XYZ = 2;
[XActual, YActual, ZActual] = InitializeXYZ(opt_XYZ);

%% Error containers initialization
Err_1  = zeros(length(rate), 6, num);
Err_2  = zeros(length(rate), 6, num);
Err_W  = zeros(length(rate), 6, num);
Err_NP = zeros(length(rate), 6, num);
Err_DK = zeros(length(rate), 6, num);

%%
for s = 1:num
    counter = 0;
    optPDF = 1;
    [A1, B1, C1, A2, B2, C2, A3, B3, C3] = ...
        generateSetsOfABC(Num, optPDF, gmean, k*Cov, XActual, YActual, ZActual);
    
    for r = rate
        counter = counter + 1;
        
        % Permutate B1, A2, C3 data streams
        [A1_perm, B1_perm, C1_perm] = permFixABC(A1(:,:,1), B1, C1, r);
        [C2_perm, A2_perm, B2_perm] = permFixABC(C2(:,:,1), A2, B2, r);
        [B3_perm, C3_perm, A3_perm] = permFixABC(B3(:,:,1), C3, A3, r);
        
        %% -- concatenate data for traditional simultaneous axbycz solvers
        A_perm = cat(3, A1_perm, A2_perm, A3_perm);
        B_perm = cat(3, B1_perm, B2_perm, B3_perm);
        C_perm = cat(3, C1_perm, C2_perm, C3_perm);
        
        %% -- different algorithms
        [X_f1,   Y_f1,   Z_f1] = axbyczProb1(A1, B1_perm, C1_perm, A2_perm, B2_perm, C2, 0, 0, 0);
        [X_f2,   Y_f2,   Z_f2] = axbyczProb2(A1, B1_perm, C1_perm, A2_perm, B2_perm, C2, A3_perm, B3, C3_perm);
        [X_wang, Y_wang, Z_wang] = Wang_AXBYCZ( A_perm, B_perm, C_perm, XActual, YActual, ZActual);
        [X_NP,   Y_NP,   Z_NP ] = Yan_AXBYCZ_PN( A_perm, B_perm, C_perm, XActual, YActual, ZActual);
        [X_DK,   Y_DK,   Z_DK ] = Yan_AXBYCZ_DK( A1, B1_perm, C1_perm, A2_perm, B2_perm, C2 );
        
        %% ----- err analysis ------
        % ------- Prob1 Error ------ %
        Err_1(counter,:,s) = getErrorAXBYCZ(X_f1, Y_f1, Z_f1, ...
            XActual, YActual, ZActual);
        % ------- Prob2 Error ------ %
        Err_2(counter,:,s) = getErrorAXBYCZ(X_f2, Y_f2, Z_f2, ...
            XActual, YActual, ZActual);
        % ------- Wang Error ------ %
        Err_W(counter,:,s) = getErrorAXBYCZ(X_wang, Y_wang, Z_wang, ...
            XActual, YActual, ZActual);
        % ------- NP Error ------ %
        Err_NP(counter,:,s) = getErrorAXBYCZ(X_NP, Y_NP, Z_NP, ...
            XActual, YActual, ZActual);
        % ------- DK Error ------ %
        Err_DK(counter,:,s) = getErrorAXBYCZ(X_DK, Y_DK, Z_DK, ...
            XActual, YActual, ZActual);
    end
end

%% Compute the averaged errors
Err1_Avg = sum(Err_1,3)/num;
Err2_Avg = sum(Err_2,3)/num;
ErrW_Avg = sum(Err_W,3)/num;
ErrNP_Avg = sum(Err_NP,3)/num;
ErrDK_Avg = sum(Err_DK,3)/num;

%% Plots
figure; hold on;
y_lb = {'$E_{\bf R_{X}}$', '$E_{\bf R_{Y}}$', '$E_{\bf R_{Z}}$',...
    '$E_{\bf t_{X}}$', '$E_{\bf t_{Y}}$', '$E_{\bf t_{Z}}$'};

for i = 1:6
    % ------- Subplot 1 ------- %
    subplot(2,3,i); hold on
    plot(rate, Err1_Avg(:,i),'d-r')
    plot(rate, Err2_Avg(:,i),'*-g')
    plot(rate, ErrW_Avg(:,i),'o-b')
    plot(rate, ErrNP_Avg(:,i),'+-m')
    plot(rate, ErrDK_Avg(:,i),'x-c')
    
    len1 = legend('Prob1','Prob2','Wang', 'PN', 'DK');
    set(len1,'FontSize',14,'Interpreter','latex','Location','northwest');
    ylabel(y_lb{i},'FontSize',18,'Interpreter','latex');
end
