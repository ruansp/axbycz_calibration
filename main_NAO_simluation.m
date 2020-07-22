% This script tests for iterative refinement on simulated data
%
% Author: Sipu Ruan, ruansp@jhu.edu, 2020

clear; close all; clc

%% Add file dependencies
addpath util/
addpath solvers/

%% Generate Data
% Initial guesses:
%  1 for identity; 2 for results from Prob 1.
init_guess = 2;
% Initial guess as Identity
if init_guess == 1
    X_init = eye(4);
    Y_init = eye(4);
    Z_init = eye(4);
end

% Choice of scramble rate
isRandPerm = 1;
r = 0:10:100;

% True values for simulations
X_true = [rotx(-pi/2)*roty(pi/2)*rotx(1.2*pi/180)...
    [58.71;0;63.64]/1e3; 0 0 0 1]; % from kinematics data
Y_true = [rotz(pi)*rotz(pi/4) [400;0;0]/1e3; 0 0 0 1];
Z_true = [rotz(pi) [0;0;10]/1e3; 0 0 0 1];

for rk = 1:size(r,2)
    length = 100;
    Mean = [0; 0; 0; 0; 0 ;0];
    Cov = 0.1*eye(6,6);
    
    A = []; B = []; C = []; Bp = [];
    
    % Generate A, B, C using different data distribution and noise
    % level
    for i = 1:5
        [A1{i}, B1{i}, C1{i}] =  generateABC(length, 1, 1, Mean, Cov, X_true, Y_true, Z_true);
        [A2{i}, B2{i}, C2{i}] =  generateABC(length, 3, 1, Mean, Cov, X_true, Y_true, Z_true);
        
        AA = cat(3, A1{i}, A2{i});
        BB = cat(3, B1{i}, B2{i});
        CC = cat(3, C1{i}, C2{i});
        
        A = cat(3, A, AA);
        B = cat(3, B, BB);
        C = cat(3, C, CC);
        
        if isRandPerm
            Bp1{i} = scrambleData(B1{i}, r(rk));
            Bp2{i} = scrambleData(B2{i}, r(rk));
        end
    end
    
    % Inputs for Prob 1
    AA1 = A1{1}; BB1 = B1{1}; CC1 = C1{1};
    AA2 = A2{1}; BB2 = B2{1}; CC2 = C2{1};
    if isRandPerm
        BBp1 = scrambleData(BB1, r(rk));
        BBp2 = scrambleData(BB2, r(rk));
        Bp = scrambleData(B, r(rk));
    end
    
    %% Prob 1
    disp('Probabilistic Method 1...')
    [X_cal1, Y_cal1, Z_cal1] = axbyczProb1(AA1(:,:,1), BBp1, CC1, ...
        AA2, BBp2, CC2(:,:,1), 0, 0, 0);
    
    % Initial guess for iterative refinement as the results from prob 1
    if init_guess == 2
        X_init = X_cal1;
        Y_init = Y_cal1;
        Z_init = Z_cal1;
    end
    
    %% Iteratice Refinement
    disp('Iteratice Refinement...')
    [X_cal2, Y_cal2, Z_cal2, num2(rk)] = axbyczProb3(A1, Bp1, C1, ...
        A2, Bp2, C2, X_init, Y_init, Z_init);
    
    %% Call traditional AXB=YCZ algorithm to solve for X, Y and Z given A, B and C
    disp('Wang Method...')
    [X_cal3, Y_cal3, Z_cal3, num3(rk)] = Wang_AXBYCZ(A, Bp, C,...
        X_init, Y_init, Z_init);
    
    %% Verification
    % Prob 1
    err_prob(rk,:) = getErrorAXBYCZ(X_cal1, Y_cal1, Z_cal1, ...
        X_true, Y_true, Z_true);
    err1(rk) = metric(A1,B1,C1,X_cal1,Y_cal1,Z_cal1) +...
        metric(A2,B2,C2,X_cal1,Y_cal1,Z_cal1);
    
    % Iterative refinement
    err_iter(rk,:) = getErrorAXBYCZ(X_cal2, Y_cal2, Z_cal2, ...
        X_true, Y_true, Z_true);
    err2(rk) = metric(A1,B1,C1,X_cal2,Y_cal2,Z_cal2) +...
        metric(A2,B2,C2,X_cal2,Y_cal2,Z_cal2);
    
    % Wang method
    err_wang(rk,:) = getErrorAXBYCZ(X_cal3, Y_cal3, Z_cal3, ...
        X_true, Y_true, Z_true);
    err3(rk) = metric(A1,B1,C1,X_cal3,Y_cal3,Z_cal3) +...
        metric(A2,B2,C2,X_cal3,Y_cal3,Z_cal3);
end

%% Plot error v.s. scramble rate
% Errors with ground truth
figure; hold on;
fontSize = 20;
lineW = 1;
y_lb = {'$R_{X}$', '$R_{Y}$', '$R_{Z}$',...
    '${\bf t_{X}}$', '${\bf t_{Y}}$', '${\bf t_{Z}}$'};

for i = 1:6
    % ------- Subplot 1 ------- %
    subplot(2,3,i); hold on
    plot(r, err_prob(:,i), 'o-r', 'LineWidth', lineW)
    plot(r, err_iter(:,i), 'd-g', 'LineWidth', lineW)
    plot(r, err_wang(:,i), '*-b', 'LineWidth', lineW)
    
    len1 = legend('Prob1','Iterative','Wang');
    ylabel(y_lb{i},'FontSize',fontSize,'Interpreter','latex');
    xlabel('Scramble Rate / %')
end

% Errors between two sides of calibration equations
figure; hold on;
fontSize = 20;
lineW = 1;

plot(r,err1, 'o-r', 'LineWidth', lineW);
plot(r,err2, 'd-g','LineWidth', lineW);
plot(r,err3, '*-b','LineWidth', lineW);

lgd = legend('Prob 1', 'Iterative', 'Wang');
lgd.FontSize = fontSize;
xlabel('Scramble Rate (%)', 'FontSize', fontSize)
ylabel('Error', 'FontSize', fontSize)

%% Metric for computing errors between two sides of calibration equations
function diff = metric(A,B,C,X,Y,Z)
diff = 0;
N = 0;
for i = 1:size(A,2)
    for j = 1:size(A{i},3)
        diff = diff + norm(A{i}(:,:,j)*X*B{i}(:,:,j)-Y*C{i}(:,:,j)*Z, 'fro');
        N = N+1;
    end
end
diff = diff/N;
end