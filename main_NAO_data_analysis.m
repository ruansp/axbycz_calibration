% This script tests for iterative refinement on real data
%
% Author: Sipu Ruan, ruansp@jhu.edu, 2017

clear; close all; clc

%% Add file dependencies
addpath util/
addpath solvers/

%% Load data
load('real_data/transform_ABC_unified.mat'); % load the experiment data
load('real_data/transform_ABC_unified_fixA.mat');
load('real_data/transform_ABC_unified_fixC.mat');

%% Generate Data
% Initial guesses:
%  1 for identity; 2(or 3) for approximate measurement from kinematics 
%  data of the robot; 3 for results from Prob 1.
init_guess = 3;
% Initial guess as Identity
if init_guess == 1
    X_init = eye(4);
    Y_init = eye(4);
    Z_init = eye(4);
else
    % Initial guess as approx. measurement
    X_init = [rotx(-pi/2)*roty(pi/2)*rotx(1.2*pi/180)...
        [58.71;0;63.64]/1e3; 0 0 0 1]; % from kinematics data
    Y_init = [rotz(pi)*rotz(pi/4) [400;0;0]/1e3; 0 0 0 1];
    Z_init = [rotz(pi) [0;0;10]/1e3; 0 0 0 1];
end

% Choice of simulation or real robot
isSim = 0;

% Choice of scramble rate
isRandPerm = 1;
r = 0:10:100;
for rk = 1:size(r,2)
    % convert cells to matrices
    A = []; B = []; C = []; Bp = [];
    for i = 1:size(headTf1,2)
        % Inputs for Iterative Refinement
        [A1{i},B1{i},C1{i}] = convertCell2Mat(headTf1{i}, handTf1{i}, tagTf1{i});
        [A2{i},B2{i},C2{i}] = convertCell2Mat(headTf2{i}, handTf2{i}, tagTf2{i});
        if isRandPerm
            Bp1{i} = scrambleData(B1{i}, r(rk));
            Bp2{i} = scrambleData(B2{i}, r(rk));
        end
        
        % Inputs for Wang
        AA = cat(3, A1{i}, A2{i});
        BB = cat(3, B1{i}, B2{i});
        CC = cat(3, C1{i}, C2{i});
        
        A = cat(3, A, AA);
        B = cat(3, B, BB);
        C = cat(3, C, CC);
    end
    
    % Inputs for Prob 1
    [AA1,BB1,CC1] = convertCell2Mat(headTf1{1}, handTf1{1}, tagTf1{1});
    [AA2,BB2,CC2] = convertCell2Mat(headTf2{1}, handTf2{1}, tagTf2{1});
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
    if init_guess == 3
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
    diff_X_rot1 = roterror(X_cal1(1:3,1:3),X_init(1:3,1:3));
    diff_Y_rot1 = roterror(Y_cal1(1:3,1:3),Y_init(1:3,1:3));
    diff_Z_rot1 = roterror(Z_cal1(1:3,1:3),Z_init(1:3,1:3));
    
    diff_X1 = norm(X_cal1(:,4)-X_init(:,4));
    diff_Y1 = norm(Y_cal1(:,4)-Y_init(:,4));
    diff_Z1 = norm(Z_cal1(:,4)-Z_init(:,4));
    
    err1(rk) = metric(A1,B1,C1,X_cal1,Y_cal1,Z_cal1) +...
        metric(A2,B2,C2,X_cal1,Y_cal1,Z_cal1);
    
    % Iterative refinement
    diff_X_rot2 = roterror(X_cal2(1:3,1:3),X_init(1:3,1:3));
    diff_Y_rot2 = roterror(Y_cal2(1:3,1:3),Y_init(1:3,1:3));
    diff_Z_rot2 = roterror(Z_cal2(1:3,1:3),Z_init(1:3,1:3));
    
    diff_X2 = norm(X_cal2(:,4)-X_init(:,4));
    diff_Y2 = norm(Y_cal2(:,4)-Y_init(:,4));
    diff_Z2 = norm(Z_cal2(:,4)-Z_init(:,4));
    
    err2(rk) = metric(A1,B1,C1,X_cal2,Y_cal2,Z_cal2) +...
        metric(A2,B2,C2,X_cal2,Y_cal2,Z_cal2);
    
    % Wang method
    diff_X_rot3 = roterror(X_cal3(1:3,1:3),X_init(1:3,1:3));
    diff_Y_rot3 = roterror(Y_cal3(1:3,1:3),Y_init(1:3,1:3));
    diff_Z_rot3 = roterror(Z_cal3(1:3,1:3),Z_init(1:3,1:3));
    
    diff_X3 = norm(X_cal3(:,4)-X_init(:,4));
    diff_Y3 = norm(Y_cal3(:,4)-Y_init(:,4));
    diff_Z3 = norm(Z_cal3(:,4)-Z_init(:,4));
    
    err3(rk) = metric(A1,B1,C1,X_cal3,Y_cal3,Z_cal3) +...
        metric(A2,B2,C2,X_cal3,Y_cal3,Z_cal3);
end

%% Plot error v.s. scramble rate
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
title('Error: Real Data', 'FontSize', fontSize)


%% Supporting functions
% Convert data to 3d matrices
function [A, B, C] = convertCell2Mat(headTf, handTf, tagTf)
for i = 1:size(headTf,2)
    A(:,:,i) = headTf{i};
    B(:,:,i) = tagTf{i};
    C(:,:,i) = handTf{i};
end
end

% Metric for computing errors
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