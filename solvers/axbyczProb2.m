function [X_final, Y_final, Z_final] = axbyczProb2(A1, B1, C1, A2, B2, C2, A3, B3, C3)
% This function implements the Prob2 method in the paper
% Prerequisites on the input
%   A1 is constant with B1 and C1 free
%   C2 is constant with A1 adn B1 free
%   B3 is constant with A3 and C3 free 
%
% Authors: Qianli Ma, qianli.ma622@gmail.com; 
%          Zachariah Goh, zach_goh@yahoo.com
% Modifications: Sipu Ruan, ruansp@jhu.edu
A1 = A1(:,:,1); C2 = C2(:,:,1); B3 = B3(:,:,1);

%% file dependencies
% addpath ../../rvctools/robot
% addpath ../../rvctools/common
% addpath ../../rvctools/screws
% addpath ../../rvctools/util
% addpath(genpath('~/Dropbox/RSS2016/Matlab/Zach/Batch'))
% addpath(genpath('~/Dropbox/RSS2016/Matlab/Zach'))

fprintf('Running Prob2 method ... \n')
%% ------ Solve for Z -------- %%
% A1 fixed, B1 and C1 free

% ------ using probability methods ------
% calculate Z_g : all guesses of Z
[ Z_g, ~, MeanC1, MeanB1, ~, ~ ] = batchSolveXY(C1, B1, false, 0, 0);

% Keep the candidates of Z that are SE3
% Normally there be four Z \in SE3
Z_index = 1;
Z = [];
for i = 1:size(Z_g,3)
    if det(Z_g(:,:,i)) > 0
        Z(:,:,Z_index) = Z_g(:,:,i);
        Z_index = Z_index + 1;
    end
end

%% ------ Solve for X -------- %%
% C2 fixed, A2 and B2 free

% ------ Calculate B2^-1 -------
Num = size(A2, 3);
A2_inv = zeros(4, 4, Num);
B2_inv = zeros(4, 4, Num);
for i = 1:Num
    A2_inv(:, :, i) = inv(A2(:, :, i));
    B2_inv(:, :, i) = inv(B2(:, :, i));
end

% ------ using probability methods ------
% calculate X_g : all guesses of X
[ X_g, ~, MeanA2, ~, ~, ~ ] = batchSolveXY(A2, B2_inv, false, 0, 0);

% Calculate MeanA2_inv and MeanB2 for the cost function later
[ ~, ~, ~, MeanB2, ~, ~ ] = batchSolveXY(A2_inv, B2, false, 0, 0);

% Keep the candidates of X that are SE3
X_index = 1;
X = [];
for i = 1:size(X_g,3)
    if det(X_g(:,:,i)) > 0
        X(:,:,X_index) = X_g(:,:,i);
        X_index = X_index + 1;
    end
end

%% ------ Solve for Y -------- %%
% B3 fixed, A3 and C3 free

% ------ Calculate B2^-1 -------
A3_inv = zeros(4, 4, Num);
C3_inv = zeros(4, 4, Num);
for i = 1:Num
    A3_inv(:, :, i) = inv(A3(:, :, i));
    C3_inv(:, :, i) = inv(C3(:, :, i));
end

% ------ using probability methods ------
% calculate X_g : all guesses of X
[ Y_g_inv, ~, ~, ~, ~, ~ ] = batchSolveXY(C3_inv, A3_inv, false, 0, 0);

% Calculate MeanA2 and MeanC2 for the cost function later
[ ~, ~, MeanC3, MeanA3, ~, ~ ] = batchSolveXY(C3, A3, false, 0, 0);

% Keep the candidates of Y that are SE3
Y_index = 1;
Y = [];
for i = 1:size(X_g,3)
    if det(Y_g_inv(:,:,i)) > 0
        Y(:,:,Y_index) = inv(Y_g_inv(:,:,i));
        Y_index = Y_index + 1;
    end
end

%% Find out the optimal (X, Y, Z) that minimizes cost
s_X = size(X, 3);
s_Y = size(Y, 3);
s_Z = size(Z, 3);

cost = zeros(s_X, s_Y*s_Z);
weight = 1.8; % weight on the translational error of the cost function

for i = 1:s_X
    for j = 1:s_Y
        for p = 1:s_Z
            
            left1  = A1*X(:,:,i)*MeanB1;
            right1 = Y(:,:,j)*MeanC1*Z(:,:,p);
            diff1  = roterror(left1, right1) + weight*tranerror(left1, right1);
            
            left2  = MeanA2*X(:,:,i)*MeanB2;
            right2 = Y(:,:,j)*C2*Z(:,:,p);
            diff2  = roterror(left2, right2) + weight*tranerror(left2, right2);
            
            left3  = MeanA3*X(:,:,i)*B3;
            right3 = Y(:,:,j)*MeanC3*Z(:,:,p);
            diff3  = roterror(left3, right3) + weight*tranerror(left3, right3);
            
            % How to better design the cost function is still an open
            % question
            cost(i, (j-1)*s_Z + p) = norm(diff1(:)) + norm(diff2(:)) + norm(diff3(:));
        end
    end
end


% recover the X,Y,Z that minimizes cost
[~, I1] = min(cost(:));

[I_row, I_col] = ind2sub(size(cost), I1);
X_final = X(:, :, I_row); % final X

if mod(I_col, s_Y) > 0
    if mode(I_col, s_Y) == s_Z
        index_Y = s_Z;
    end
    index_Y = floor(I_col/s_Y) + 1;
else
    index_Y = floor(I_col/s_Y);
end

Y_final = Y(:, :, index_Y); % final Y

if rem(I_col,4) > 0
    index_Z = rem(I_col, 4) ;
else
    index_Z = 4;
end

Z_final = Z(:, :, index_Z); % final Z

end