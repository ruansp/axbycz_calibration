function [ X, Y, Z ] = Yan_AXBYCZ_DK( A1, B1, C1, A2, B2, C2 )
% Implements the D-K algorithm in Yan et al (2015)
%   Solves for X, Y, Z in the matrix equation AXB=YCZ given A, B, C

% Input: A1 and C2 are 4x4 homogeneous matrices since they are fixed 
%        B1,C1, A2,B2 are 4 x 4 x n homogeneous matrices
% Output: X, Y, Z are 4 x 4 homogeneous matrices

A1 = A1(:,:,1); C2 = C2(:,:,1);

num = size(B1, 3); % number of measurements
fprintf('Num of data: %d\n', num); 

%% Solve AX=YB type of equations using Li(2010) method of Kronecker product
[Z, Xt] = Li_AXYB_kron( C1, B1 ); % fixing A
% [Z, Xt] = li(C1,B1); % use shah's implementation
Binv = zeros(size(B2)); % 4x4xnum
for i=1:num
  Binv(:,:,i) = inv(B2(:,:,i)); % calc the inv of B b4 passing to solver
end
[X, Zt] = Li_AXYB_kron( A2, Binv ); % fixing C
% [X, Zt] = li(A2,Binv); % use shah's implementation

%% Enforce orthogonality on rotational part of Xt and Z, X and Zt
%  by finding the nearest orthogonal matrix using SVD
[U,~,V] = svd(Xt(1:3,1:3));
Xt(1:3,1:3) = U*V';
[U,~,V] = svd(X(1:3,1:3));
X(1:3,1:3) = U*V';
[U,~,V] = svd(Zt(1:3,1:3));
Zt(1:3,1:3) = U*V';
[U,~,V] = svd(Z(1:3,1:3));
Z(1:3,1:3) = U*V';

%% calc the two possibilities of Y and choose the one with the smallest
% error
Y1 = A1 * X / Xt;
Y2 = Zt / Z / C2;

Err1 = norm( A1*X*B1(:,:,1) - Y1*C1(:,:,i)*Z, 'fro' );
Err2 = norm( A2(:,:,1)*X*B2(:,:,1) - Y2*C2*Z, 'fro' );
if ( Err1 < Err2 )
  Y = Y1;
else
  Y = Y2;
end
end

%% implements the Kronecker product method in Li et al (2010) paper
function [X, Y] = Li_AXYB_kron( A, B )

num = size(A, 3); % number of measurements

RA = A(1:3, 1:3, :); % 3x3xnum
RB = B(1:3, 1:3, :);
TA = A(1:3, 4, :);  % 3x1xnum
TB = B(1:3, 4, :);

K = zeros(12*num, 24);
t = zeros(12*num, 1);
for i = 1:num
  K( (i-1)*12+1:(i-1)*12+9, 1:9 ) = kron( RA(:,:,i), eye(3) );
  K( (i-1)*12+1:(i-1)*12+9, 10:18 ) = -kron( eye(3), RB(:,:,i)' );
  K( (i-1)*12+10:i*12, 10:18 ) = kron(eye(3), TB(:,:,i)' );
  K( (i-1)*12+10:i*12, 19:21 ) = -RA(:,:,i);
  K( (i-1)*12+10:i*12, 22:24 ) = eye(3);
  t( (i-1)*12+10:i*12 ) = TA(:,:,i);
end

% solve Kv = t using least squares
v = pinv(K) * t; % 24x1

X = zeros(4,4); X(4,4) = 1;
Y = zeros(4,4); Y(4,4) = 1;

% reform the X, Y homogeneous matrices from the vectorized versions
X(1:3,1:3) = reshape(v(1:9), 3, 3)'; %need transpose because reshape goes down col first then row
Y(1:3,1:3) = reshape(v(10:18), 3, 3)';
X(1:3, 4) = reshape(v(19:21), 3, 1);
Y(1:3, 4) = reshape(v(22:24), 3, 1);
end


