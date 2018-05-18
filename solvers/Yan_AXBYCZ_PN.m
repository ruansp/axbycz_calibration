function [ X, Y, Z ] = Yan_AXBYCZ_PN( A, B, C, Xact, Yact, Zact)
% Implements the D-K algorithm in Yan et al (2015)
%   Solves for X, Y, Z in the matrix equation AXB=YCZ given A, B, C

% Input: A1 and C2 are 4x4 homogeneous matrices since they are fixed
%        B1,C1, A2,B2 are 4 x 4 x n homogeneous matrices
% Output: X, Y, Z are 4 x 4 homogeneous matrices

num = size(A, 3); % number of measurements
fprintf('Num of data: %d\n', num);

%% Get rotation and translation components of A, B, C
RA = A(1:3, 1:3, :); % 3x3xnum
RB = B(1:3, 1:3, :);
RC = C(1:3, 1:3, :);
TA = A(1:3, 4, :);  % 3x1xnum
TB = B(1:3, 4, :);
TC = C(1:3, 4, :);

%% Form error function - nested because it needs to access some variables above
    function F = myfun(xyz)
        % xyz is a 21x1 vector of [qx qy qz tx ty tz]
        F = zeros(num*12,1);
        
        qx = Quaternion(xyz(1:4)).unit(); % normalize quaternion
        RX = qx.R;
        qy = Quaternion(xyz(5:8)).unit();
        RY = qy.R;
        qz = Quaternion(xyz(9:12)).unit();
        RZ = qz.R;
        TX = xyz(13:15);
        TY = xyz(16:18);
        TZ = xyz(19:21);
        
        for i = 1:num
            Rerr = RA(:,:,i)*RX*RB(:,:,i)-RY*RC(:,:,i)*RZ; % 3x3
            Terr = RA(:,:,i)*RX*TB(:,:,i)+RA(:,:,i)*TX+TA(:,:,i)...
                - RY*RC(:,:,i)*TZ-RY*TC(:,:,i)-TY; % 3x1
            F( (i-1)*12+1 : (i-1)*12+9 ) = reshape(Rerr, 9, 1);
            F( (i-1)*12+10 : (i-1)*12+12 ) = Terr;
        end
        
    end

%% get initial estimate of quaternions and translation
opt  = 2;

if opt == 1
    % a) Use perturbation of actual X, Y, Z
    e = pi/5;
    RX_init = expm( so3_vec(e*ones(3,1)) ) * Xact(1:3,1:3); % 3x3
    RY_init = expm( so3_vec(e*ones(3,1)) ) * Yact(1:3,1:3);
    RZ_init = expm( so3_vec(e*ones(3,1)) ) * Zact(1:3, 1:3);
    TX_init = Xact(1:3, 4) + e*ones(3,1);  % 3x1xnum
    TY_init = Yact(1:3, 4) + e*ones(3,1);
    TZ_init = Zact(1:3, 4) + e*ones(3,1);
    
elseif opt == 2
    % b) Use randomly generated X, Y, Z
    M = zeros(6,1); %mean
    Sig = eye(6)*2; %covariance
    XActual = expm(se3_vec(mvg(M, Sig, 1)));
    YActual = expm(se3_vec(mvg(M, Sig, 1)));
    ZActual = expm(se3_vec(mvg(M, Sig, 1)));
    RX_init = XActual(1:3,1:3); % 3x3
    RY_init = YActual(1:3,1:3);
    RZ_init = ZActual(1:3,1:3);
    TX_init = XActual(1:3,4);  % 3x1xnum
    TY_init = YActual(1:3,4);
    TZ_init = ZActual(1:3,4);
    %--------------------------
end
qx_init = rotm2quat(RX_init);
qy_init = rotm2quat(RX_init);
qz_init = rotm2quat(RX_init);

% qx_init = Quaternion(RX_init).double();
% qy_init = Quaternion(RY_init).double();
% qz_init = Quaternion(RZ_init).double();
x0 = [qx_init'; qy_init'; qz_init'; TX_init; TY_init; TZ_init];

%% call lsqnonlin() using Levenberg-Marquardt algorithm
options = optimoptions(@lsqnonlin,'Algorithm','levenberg-marquardt');
[res,resnorm,~,~,output] = lsqnonlin(@myfun,x0,[],[],options);

%% reform X, Y, Z from "res"
X = zeros(4); Y = zeros(4); Z = zeros(4);
X(4,4) = 1;   Y(4,4) = 1;   Z(4,4) = 1;

X(1:3,1:3) = Quaternion(res(1:4)).unit().R; %normalize b4 converting to rot matrix
Y(1:3,1:3) = Quaternion(res(5:8)).unit().R;
Z(1:3,1:3) = Quaternion(res(9:12)).unit().R;
X(1:3,4) = res(13:15); % translation components
Y(1:3,4) = res(16:18);
Z(1:3,4) = res(19:21);
end

