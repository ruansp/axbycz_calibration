function [ X, Y, Z, n_step ] = Wang_AXBYCZ( A, B, C, Xact, Yact, Zact)
% Implements the algorithm in Wang et al (2014)
%   Solves for X, Y, Z in the matrix equation AXB=YCZ given A, B, C

% Input: A, B, C are 4 x 4 x n homogeneous matrices
%       Xact, Yact, Zact are used to get a good initial estimate
% Output: X, Y, Z are 4 x 4 homogeneous matrices

num = size(A, 3); % number of measurements
% fprintf('Num of data: %d\n', num); 

%% Get rotation and translation components of A, B, C
RA = A(1:3, 1:3, :); % 3x3xnum
RB = B(1:3, 1:3, :);
RC = C(1:3, 1:3, :);
TA = A(1:3, 4, :);  % 3x1xnum
TB = B(1:3, 4, :);
TC = C(1:3, 4, :);


%% ============ Solve for RX, RY, RZ first ==============
% Set initial guess of RX, RY, RZ
e = pi/5;
RX_init = expm( so3_vec(e*ones(3,1)) ) * Xact(1:3,1:3);
RZ_init = expm( so3_vec(e*ones(3,1)) ) * Zact(1:3,1:3);
RY_init = RA(:,:,1) * RX_init * RB(:,:,1) / RZ_init / RC(:,:,1);

% figure
% trplot(RY_init, 'color', 'b');
%     hold on
%     trplot(Yact(1:3,1:3), 'color', 'r'); 
    
% fprintf('Err in initial guess: RX = %.5f, RY = %.5f, RZ = %.5f\n',...
%   roterror( RX_init, Xact(1:3,1:3) ), ...
%   roterror( RY_init, Yact(1:3,1:3) ),...
%   roterror( RZ_init, Zact(1:3,1:3) ) );

% Iterate until norm of delR = [delRX; delRY; delRZ] falls below a predefined threshold
delR = 10000 * ones(9,1);  % use a large value initially 

n_step = 0;

while (norm(delR) > .01 && n_step < 500 )
  
  q = zeros(num*9,1); % q_tilde in paper
  F = zeros(num*9,9); % F_tilde in paper
  
  for i = 1:num
    tmp1 = RX_init * RB(:,:,i);
    tmp2 = RY_init * RC(:,:,i) * RZ_init;
    qq = -RA(:,:,i) * tmp1 + tmp2;
    q( (i-1)*9+1:i*9 ) = [qq(:,1); qq(:,2); qq(:,3)]; % 9x1
    
    F11 = -RA(:,:,i) * so3_vec( tmp1(:,1) ); % 3x3
    F21 = -RA(:,:,i) * so3_vec( tmp1(:,2) );
    F31 = -RA(:,:,i) * so3_vec( tmp1(:,3) );
    F12 = so3_vec( tmp2(:,1) );
    F22 = so3_vec( tmp2(:,2) );
    F32 = so3_vec( tmp2(:,3) );
    F13 = RY_init * RC(:,:,i) * so3_vec( RZ_init(:,1) );
    F23 = RY_init * RC(:,:,i) * so3_vec( RZ_init(:,2) );
    F33 = RY_init * RC(:,:,i) * so3_vec( RZ_init(:,3) );
    F( (i-1)*9+1:i*9, : ) = [ F11 F12 F13;
                              F21 F22 F23;
                              F31 F32 F33 ]; 
  end
  
  delR = (F'*F) \ F' * q;  % = inv(F'F)F'q
  
%   fprintf('cost is %d \n', norm(delR))
  
  thetaX = norm( delR(1:3) );
  RX_init = skewexp( delR(1:3)/thetaX, thetaX ) * RX_init;
  
  thetaY = norm( delR(4:6) );
  RY_init = skewexp( delR(4:6)/thetaY, thetaY ) * RY_init;
  
  thetaZ = norm( delR(7:9) );
  RZ_init = skewexp( delR(7:9)/thetaZ, thetaZ ) * RZ_init;
  
  n_step = n_step + 1;
end

%% ============ Solve for TX, TY, TZ next ==============
J = zeros(3*num, 9); % J_tilde
p = zeros(3*num, 1); % p_tilde

for i = 1:num
  J( ((i-1)*3+1):(i*3), : ) = [ RA(:,:,i)  -eye(3)  -RY_init * RC(:,:,i) ];
  p( ((i-1)*3+1):(i*3) ) = -TA(:,:,i) - RA(:,:,i)*RX_init*TB(:,:,i) + RY_init*TC(:,:,i);
end

translation = (J'*J) \ J' * p;
tX = translation(1:3);
tY = translation(4:6);
tZ = translation(7:9);

%% Form the homogeneous matrices for X, Y, Z
X = zeros(4); Y = zeros(4); Z = zeros(4);
X(4,4) = 1;   Y(4,4) = 1;   Z(4,4) = 1;
X(1:3,1:3) = RX_init;
Y(1:3,1:3) = RY_init;
Z(1:3,1:3) = RZ_init;
X(1:3,4) = tX;
Y(1:3,4) = tY;
Z(1:3,4) = tZ;
end

