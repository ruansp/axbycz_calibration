function [X_final, Y_final, Z_final] = axbyczProb1(A1, B1, C1, A2, B2, C2, opt, nstd1, nstd2)
% This function implements the Prob1 method in the paper
% Prerequisites on the input
%   A1 is constant with B1 and C1 free
%   C2 is constant with A2 adn B2 free
%
% Authors: Qianli Ma, qianli.ma622@gmail.com; 
%          Zachariah Goh, zach_goh@yahoo.com
% Modifications: Sipu Ruan, ruansp@jhu.edu

A1 = A1(:,:,1); C2 = C2(:,:,1);

% fprintf('Running Prob1 method ... \n')
%% ------ Solve for Z -------- %%
% A1 fixed, B1 and C1 free

% ------ using probability methods ------
% calculate Z_g : all guesses of Z
[ Z_g, ~, MeanC1, MeanB1, ~, ~ ] = batchSolveXY(C1, B1, opt, nstd1, nstd2);

% Keep the candidates of Z that are SE3
% Normally, there will be four Z \in SE3
Z_index = 1;
Z = [];
for i = 1:size(Z_g,3)
    if det(Z_g(:,:,i)) > 0
        Z(:,:,Z_index) = Z_g(:,:,i);
        Z_index = Z_index + 1;
    end
end

s_Z = size(Z,3);

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
[ X_g, ~, MeanA2, ~, ~, ~ ] = batchSolveXY(A2, B2_inv, opt, nstd1, nstd2);

% Calculate MeanB2 for computing Y later
% Note: can be further simplied by using only the distribution function
[ ~, ~, ~, MeanB2, ~, ~ ] = batchSolveXY(A2_inv, B2, opt, nstd1, nstd2);

% Keep the candidates of X that are SE3
% Normally, there will be four X \in SE3
X_index = 1;
X = [];
for i = 1:size(X_g,3)
    if det(X_g(:,:,i)) > 0
        X(:,:,X_index) = X_g(:,:,i);
        X_index = X_index + 1;
    end
end

s_X = size(X, 3);

%% ------ Solve for Y -------- %%
% Compute Y using the mean equations
Y = zeros(4, 4, 2*s_X*s_Z);
for i = 1:s_X
    for j = 1:s_Z
        
        % There are at least four mean equations to choose from to compute
        % Y. It will be interesting to see how each choice of the mean 
        % equations can affect the result
        Y(:,:,(i-1)*s_Z + j) = (A1.*X(:,:,i)*MeanB1/Z(:,:,j)) / MeanC1;
        Y(:,:,(i-1)*s_Z + j +s_X*s_Z) = (MeanA2*X(:,:,i)*MeanB2 / Z(:,:,j)) / C2;
    end
end

%% Find out the optimal (X, Y, Z) that minimizes cost

s_Y = size(Y, 3);

cost = zeros(s_X, s_Y*s_Z);
weight = 1.5; % weight on the translational error of the cost function
for i = 1:s_X
    for j = 1:s_Z
        for m = 1:s_Y
            
            left1  = A1*X(:,:,i)*MeanB1;
            right1 = Y(:,:,m)*MeanC1*Z(:,:,j);

            diff1 = roterror(left1, right1) + weight*tranerror(left1, right1);
            
            left2  = MeanA2*X(:,:,i)*MeanB2;
            right2 = Y(:,:,m)*C2*Z(:,:,j);
            diff2 = roterror(left2, right2) + weight*tranerror(left2, right2);
            
            % different error metrics can be picked and this (diff1 + 
            % diff2) is the best one so far. However, it can still be 
            % unstable sometimes, and miss the optimal solutions            
            cost(i,(j-1)*s_Y+m) = norm(diff1(:)) + norm(diff2(:));
        end
    end
end

% recover the X,Y,Z that minimizes cost1
[~, I1] = min(cost(:));


[I_row, I_col] = ind2sub(size(cost), I1);
X_final = X(:, :, I_row); % final X

if mod(I_col, s_Y) > 0
    index_Z = floor(I_col/s_Y) + 1;
else
    index_Z = floor(I_col/s_Y);
end

Z_final = Z(:, :, index_Z); % final Z


if rem(I_col,s_Y) > 0
    index_Y = rem(I_col, s_Y) ;
else
    index_Y = s_Y;
end

Y_final = Y(:, :, index_Y); % final Y


end