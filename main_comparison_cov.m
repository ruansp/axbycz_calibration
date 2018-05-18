% This is the main file for probabilistic AXB = YCZ problem
%  Comparisons with covariance of noise
%
% Authors: Qianli Ma, qianli.ma622@gmail.com; 
%          Zachariah Goh, zach_goh@yahoo.com
% Modifications: Sipu Ruan, ruansp@jhu.edu

clear; close all; clc

%% Add file dependencies
addpath util/
addpath solvers/

%% Initialize Parameters
counter = 0;
gmean = [0; 0; 0; 0; 0 ;0];
coeff1 = 0.02:0.02:0.14; % Scaling factor for the covariances
cov = eye(6,6); % coeff1*cov
coeff2 = 0.02:0.04:0.14; % Scaling factor for the covariances
cov_noise = eye(6,6); % coeff2*cov
num = 10; % number of simulations
Num = 100; % number of data
optPlot = 'lineplot'; % Plot the averaged error : 'lineplot' & ''boxplot'

opt_XYZ = 2; % generate random X, Y and Z
[XActual, YActual, ZActual] = InitializeXYZ(opt_XYZ);

% Error container initialization
if strcmp(optPlot, 'boxplot')
    Err11  = zeros(length(coeff1), num, 6);
    Err21  = zeros(length(coeff1), num, 6);    
elseif strcmp(optPlot, 'lineplot')
    Err11  = zeros(length(coeff1), 6, num);
    Err21  = zeros(length(coeff1), 6, num);
    Err12  = zeros(length(coeff1), 6, num);
    Err22  = zeros(length(coeff1), 6, num);
    Err13  = zeros(length(coeff1), 6, num);
    Err23  = zeros(length(coeff1), 6, num);
end
%%
for k = coeff1
    
    counter = counter + 1;
    
    for s = 1:num
        %% Generate data triples with different distributions
        optPDF = 1;
        [A11, B11, C11, A21, B21, C21, A31, B31, C31] = ...
            generateSetsOfABC(Num, optPDF, gmean, k*cov, XActual, YActual, ZActual);
        
        optPDF = 2;
        [A12, B12, C12, A22, B22, C22, A32, B32, C32] = ...
            generateSetsOfABC(Num, optPDF, gmean, k*cov, XActual, YActual, ZActual);
        
        optPDF = 3;
        [A13, B13, C13, A23, B23, C23, A33, B33, C33] = ...
            generateSetsOfABC(Num, optPDF, gmean, k*cov, XActual, YActual, ZActual);
        
        %% Solve for X, Y and Z
        [X_f11, Y_f11, Z_f11] = axbyczProb1(A11, B11, C11, A21, B21, C21, 0, 0, 0);
        [X_f21, Y_f21, Z_f21] = axbyczProb2(A11, B11, C11, A21, B21, C21, ...
                                            A31, B31, C31);
        
        [X_f12, Y_f12, Z_f12] = axbyczProb1(A12, B12, C12, A22, B22, C22, 0, 0, 0);
        [X_f22, Y_f22, Z_f22] = axbyczProb2(A12, B12, C12, A22, B22, C22, ...
                                            A32, B32, C32);
                                        
        [X_f13, Y_f13, Z_f13] = axbyczProb1(A13, B13, C13, A23, B23, C23, 0, 0, 0);
        [X_f23, Y_f23, Z_f23] = axbyczProb2(A13, B13, C13, A23, B23, C23, ...
                                            A33, B33, C33);
        
        %% ----- Error Analysis ------
        if strcmp(optPlot, 'boxplot')
            % ------- Prob1 Error ------ %
            Err11(counter,s,:) = getErrorAXBYCZ(X_f1, Y_f1, Z_f1, ...
                                                XActual, YActual, ZActual);
        elseif strcmp(optPlot, 'lineplot')
            
            % ------- Prob1 Error with Data of 1st Distribution ------ %
            Err11(counter,:,s) = getErrorAXBYCZ(X_f11, Y_f11, Z_f11, ...
                                                XActual, YActual, ZActual);
                                            
            % ------- Prob1 Error with Data of 2nd Distribution ------ %
            Err12(counter,:,s) = getErrorAXBYCZ(X_f12, Y_f12, Z_f12, ...
                                                XActual, YActual, ZActual);

            % ------- Prob1 Error with Data of 3rd Distribution ------ %
            Err13(counter,:,s) = getErrorAXBYCZ(X_f13, Y_f13, Z_f13, ...
                                                XActual, YActual, ZActual);
                                            
            % ------- Prob2 Error with Data of 1st Distribution ------ %
            Err21(counter,:,s) = getErrorAXBYCZ(X_f21, Y_f21, Z_f21, ...
                                                XActual, YActual, ZActual);
                                            
            % ------- Prob2 Error with Data of 2nd Distribution ------ %
            Err22(counter,:,s) = getErrorAXBYCZ(X_f22, Y_f22, Z_f22, ...
                                                XActual, YActual, ZActual);
                                            
            % ------- Prob2 Error with Data of 3rd Distribution ------ %
            Err23(counter,:,s) = getErrorAXBYCZ(X_f23, Y_f23, Z_f23, ...
                                                XActual, YActual, ZActual);
                                            
        end
        
        %% ------ plot X and Y ------
        if false
            try
                figure((counter-1)*3 + 1);
                trplot(XActual(:,:), 'color', 'r');
                hold on
                trplot(X_f11(:,:),'color','b');
                axis auto
                hold on
                title('Final X')

                figure((counter-1)*3 + 2);
                trplot(YActual(:,:), 'color', 'r');
                hold on
                trplot(Y_f11(:,:),'color','b');
                axis auto
                hold on
                title('Final Y')

                figure((counter-1)*3 + 3);
                trplot(ZActual(:,:), 'color', 'r');
                hold on
                trplot(Z_f11(:,:),'color','b');
                axis auto
                hold on
                title('Final Z')

            catch
                display(Y)
            end
        end
    end
    
end

%% Plot the averaged error w.r.t. covariances
plotProbResults(Err11, Err21, coeff1, optPlot)
plotProbResults(Err12, Err22, coeff1, optPlot)
plotProbResults(Err13, Err23, coeff1, optPlot)

