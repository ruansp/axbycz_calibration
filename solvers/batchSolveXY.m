function [ X, Y, MeanA, MeanB, SigA, SigB ] = batchSolveXY(A, B, opt, nstd_A, nstd_B)
%% Calculate X using mean and covariance propogation
% calculate Y using mean propogation (AX = YB)
%
% Input: A : 4 x 4 x n
%        B : 4 x 4 x n
%        opt : true - known nstd of data noise 
%        nstd_A : standard deviation of A
%        nstd_B : standard deviation of B
%           
% Output:
%       X : 4 x 4 x 8
%       Y : 4 x 4 x 8

%%
X_candidate = zeros(4,4,8);
Y_candidate = zeros(4,4,8);

% Reshape A and B for matching the input sizes of mex functions
[a1, a2, a3] = size(A);
A_mex = reshape(A, a1, a2*a3);
B_mex = reshape(B, a1, a2*a3);

% non-mex function version
% [ MeanA, SigA ] = distibutionPropsMex(A_mex);
% [ MeanB, SigB ] = distibutionPropsMex(B_mex);
[MeanA, SigA] = meanCov(A);
[MeanB, SigB] = meanCov(B);

% update SigA and SigB if nstd_A and nstd_B are known
if opt
    SigA = SigA - nstd_A*eye(6, 6);
    SigB = SigB - nstd_B*eye(6, 6);
end

[ VA, ~ ] = eig( SigA(1:3,1:3) );
[ VB, ~ ] = eig( SigB(1:3,1:3) );

Q1 = eye(3);
Q2 = [-1 0 0; 0 -1 0; 0 0  1];
Q3 = [-1 0 0; 0  1 0; 0 0 -1];
Q4 = [ 1 0 0; 0 -1 0; 0 0 -1];

Rx_solved = zeros(3,3,8);

% There are 8 possiblities of Rx
Rx_solved(:,:,1) = VA*Q1*VB';
Rx_solved(:,:,2) = VA*Q2*VB';
Rx_solved(:,:,3) = VA*Q3*VB';
Rx_solved(:,:,4) = VA*Q4*VB';
Rx_solved(:,:,5) = VA*-Q1*VB';
Rx_solved(:,:,6) = VA*-Q2*VB';
Rx_solved(:,:,7) = VA*-Q3*VB';
Rx_solved(:,:,8) = VA*-Q4*VB';

for i = 1:1:8
    tx_temp = ... 
        so3_vec(((Rx_solved(:,:,i)'*SigA(1:3,1:3)*Rx_solved(:,:,i))^-1*...
        (SigB(1:3,4:6)-Rx_solved(:,:,i)'*SigA(1:3,4:6)*Rx_solved(:,:,i)))');
    tx = -Rx_solved(:,:,i)*tx_temp;
    X_candidate(:,:,i) = [Rx_solved(:,:,i) tx; [0 0 0] 1];
    Y_candidate(:,:,i) = MeanA * X_candidate(:,:,i) /(MeanB);
end

X = X_candidate;
Y = Y_candidate;

end