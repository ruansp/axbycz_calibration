function [X_cal, Y_cal, Z_cal, num] = axbyczProb3(A1, B1, C1, A2, B2, C2, Xinit, Yinit, Zinit)
% This function provides iterative refinement for probabilistic solvers for
%  the AXB=YCZ calibration problem.
% It solves the system of equations:
%  (1)  A_i X \bar{B_i} = Y \bar{B_i} Z
%  (2)  \Sigma^1_{B_i} = R_Z^T \Sigma^1_{C_i} R_Z
%  (3)  C_j Z \bar{B_j}^{-1} = Y^{-1} \bar{A_j} X
%  (4)  R_{\bar{B_j}}^T \Sigma^1_{B_j} R_{\bar{B_j}}^T = R_X^T \Sigma^1_{A_j} R_X
% where X = X_init(I+\xi_X), Y = Y_init(I+\xi_Y), Z = Z_init(I+\xi_Z)
%
% Inputs:
%   A1,B1,C1: Cell arrays, which stores when fixing A1 at different poses;
%   A2,B2,C2: Cell arrays, which stores when fixing C2 at different poses;
%   Xinit,Yinit,Zinit: Initial guesses of X,Y,Z matrices.
% Outputs:
%   X_cal,Y_cal,Z_cal: Calibrated X,Y,Z matrices
%
% Updates:
%   Used the whole covariance matrix
%
% Author: Sipu Ruan, ruansp@jhu.edu, November 2017

%% Initiation
Ni = size(A1,2);
Nj = size(C2,2);
X_cal = Xinit; Y_cal = Yinit; Z_cal = Zinit;

Xupdate = Xinit; Yupdate = Yinit; Zupdate = Zinit;
xi = ones(18,1);
max_num = 500;
tol = 1e-5;
num = 1;

%% Calculate mean and covariance of varying data
for i = 1:Ni
    [A1_m(:,:,i), SigA1(:,:,i)] = meanCov(A1{i});
    [B1_m(:,:,i), SigB1(:,:,i)] = meanCov(B1{i});
    [C1_m(:,:,i), SigC1(:,:,i)] = meanCov(C1{i});
end

% invert B2
for i = 1:size(B2,2)
    for j = 1:size(B2{i},3)
        B2inv{i}(:,:,j) = inv(B2{i}(:,:,j));
    end
end

for j = 1:Nj
    [A2_m(:,:,j), SigA2(:,:,j)] = meanCov(A2{j});
    [B2_m(:,:,j), SigB2(:,:,j)] = meanCov(B2{j});
    [B2inv_m(:,:,j), SigB2inv(:,:,j)] = meanCov(B2inv{j});
    [C2_m(:,:,j), SigC2(:,:,j)] = meanCov(C2{j});
end

%% Calculate M and b matrices when fixing A and C separately
diff = metric(A1,B1,C1,Xupdate,Yupdate,Zupdate) + ...
    metric(A2,B2,C2,Xupdate,Yupdate,Zupdate);
diff = 1;
while norm(xi) >= tol && diff >= tol && num <= max_num
    for i = 1:Ni
        [MM{i}, bb{i}] = MbMat_1(A1_m(:,:,i), Xupdate, B1_m(:,:,i),...
            Yupdate, C1_m(:,:,i), Zupdate,...
            SigB1(:,:,i), SigC1(:,:,i));
    end
    
    for j = 1:Nj
        [MM{j+Ni}, bb{j+Ni}] = MbMat_2(C2_m(:,:,j), Zupdate, B2inv_m(:,:,j),...
            SE3inv(Yupdate), A2_m(:,:,j), Xupdate,...
            SigB2(:,:,j), SigA2(:,:,j), B2_m(:,:,j));
    end
    
    % Concatenate M and b matrices together
    M = []; b = [];
    M1 = [];
    M2 = [];
    M3 = [];
    M4 = [];
    for k = 1:Ni+Nj
        M = [M; MM{k}];
        b = [b; bb{k}];
    end
    for k = 1:Ni
        M1 = [M1; MM{k}(1:12,:)];
        M2 = [M2; MM{k}(13:21,:)];
    end
    for k = Ni+1:Ni+Nj
        M3 = [M3; MM{k}(1:12,:)];
        M4 = [M4; MM{k}(13:21,:)];
    end
%     rank(M)
%     rank(M1)
%     rank(M2)
%     rank(M3)
%     rank(M4)
%     rank([M1;M2])
    
    %% Inversion to get \xi_X, \xi_Y, \xi_Z
    xi = (M'*M) \ (M'*b);
%     xi = pinv(M) * b;
%     xi = lsqminnorm(M,b);

%     disp(['cost is: ', num2str(norm(xi))]);
    
    for i = 1:Ni
        diff1 = norm(A1_m(:,:,i)*Xupdate*B1_m(:,:,i)-Yupdate*C1_m(:,:,i)*Zupdate,'fro');
    end
    for i = 1:Nj
        diff2 = norm(A2_m(:,:,i)*Xupdate*B2_m(:,:,i)-Yupdate*C2_m(:,:,i)*Zupdate,'fro');
    end
    
    w_X = xi(1:3); v_X = xi(4:6);
    w_Y = xi(7:9); v_Y = xi(10:12);
    w_Z = xi(13:15); v_Z = xi(16:18);
    
    X_hat = [skew(w_X), v_X; zeros(1,4)];
    Y_hat = [skew(w_Y), v_Y; zeros(1,4)];
    Z_hat = [skew(w_Z), v_Z; zeros(1,4)];
    
%     X_cal = Xupdate*(eye(4)+X_hat);
%     Y_cal = Yupdate*(eye(4)+Y_hat);
%     Z_cal = Zupdate*(eye(4)+Z_hat);
    
    X_cal = Xupdate*expm(X_hat);
    Y_cal = Yupdate*expm(Y_hat);
    Z_cal = Zupdate*expm(Z_hat);
    
    %% Update
    Xupdate = X_cal;
    Yupdate = Y_cal;
    Zupdate = Z_cal;
    num = num+1;
    
    %% Error
    diff = metric(A1,B1,C1,Xupdate,Yupdate,Zupdate) + ...
           metric(A2,B2,C2,Xupdate,Yupdate,Zupdate);
%     disp(['diff is: ', num2str(diff)]);
end
% disp(['Number of iterations: ', num2str(num)])
end


function [M, b] = MbMat_1(A,X,B,Y,C,Z,SigB,SigC)
% Construction M and b matrices
e1 = [1;0;0]; e2 = [0;1;0]; e3 = [0;0;1];

%% AXB=YCZ
% Rotation part
M11 = -A(1:3,1:3) * X(1:3,1:3) * skew(B(1:3,1:3)*e1);
M13 = Y(1:3,1:3) * skew(C(1:3,1:3)*Z(1:3,1:3)*e1);
M15 = Y(1:3,1:3) * C(1:3,1:3) * Z(1:3,1:3) * skew(e1);

M21 = -A(1:3,1:3) * X(1:3,1:3) * skew(B(1:3,1:3)*e2);
M23 = Y(1:3,1:3) * skew(C(1:3,1:3)*Z(1:3,1:3)*e2);
M25 = Y(1:3,1:3) * C(1:3,1:3) * Z(1:3,1:3) * skew(e2);

M31 = -A(1:3,1:3) * X(1:3,1:3) * skew(B(1:3,1:3)*e3);
M33 = Y(1:3,1:3) * skew(C(1:3,1:3)*Z(1:3,1:3)*e3);
M35 = Y(1:3,1:3) * C(1:3,1:3) * Z(1:3,1:3) * skew(e3);

% Translation part
M41 = -A(1:3,1:3) * X(1:3,1:3) * skew(B(1:3,4));
M42 = A(1:3,1:3) * X(1:3,1:3);
M43 = Y(1:3,1:3) * skew(C(1:3,1:3)*Z(1:3,4) + C(1:3,4));
M44 = -Y(1:3,1:3);
M46 = -Y(1:3,1:3) * C(1:3,1:3) * Z(1:3,1:3);

M = [M11, zeros(3,3), M13, zeros(3,3),        M15, zeros(3,3);
    M21, zeros(3,3), M23, zeros(3,3),        M25, zeros(3,3);
    M31, zeros(3,3), M33, zeros(3,3),        M35, zeros(3,3);
    M41,        M42, M43,        M44, zeros(3,3),       M46];

% RHS
RHS = - A * X * B + Y * C * Z;

b = [RHS(1:3,1); RHS(1:3,2); RHS(1:3,3); RHS(1:3,4)];

%% SigBi = Ad^{-1}(Z) * SigCi * Ad^{-T}(Z)
% First block
M55 = -skew(SigB(1:3,1)) + SigB(1:3,1:3) * skew(e1);
M56 = zeros(3,3);
M65 = -skew(SigB(1:3,2)) + SigB(1:3,1:3) * skew(e2);
M66 = zeros(3,3);
M75 = -skew(SigB(1:3,3)) + SigB(1:3,1:3) * skew(e3);
M76 = zeros(3,3);

% Second block
M85 = -skew(SigB(1:3,4)) + SigB(1:3,4:6) * skew(e1);
M86 = SigB(1:3,1:3) * skew(e1);
M95 = -skew(SigB(1:3,5)) + SigB(1:3,4:6) * skew(e2);
M96 = SigB(1:3,1:3) * skew(e2);
M105 = -skew(SigB(1:3,6)) + SigB(1:3,4:6) * skew(e3);
M106 = SigB(1:3,1:3) * skew(e3);

% Third block
M115 = -skew(SigB(4:6,1)) + SigB(4:6,1:3) * skew(e1);
M116 = -skew(SigB(1:3,1));
M125 = -skew(SigB(4:6,2)) + SigB(4:6,1:3) * skew(e2);
M126 = -skew(SigB(1:3,2));
M135 = -skew(SigB(4:6,3)) + SigB(4:6,1:3) * skew(e3);
M136 = -skew(SigB(1:3,3));

% Fourth block
M145 = -skew(SigB(4:6,4)) + SigB(4:6,4:6) * skew(e1);
M146 = -skew(SigB(1:3,4)) + SigB(4:6,1:3) * skew(e1);
M155 = -skew(SigB(4:6,5)) + SigB(4:6,4:6) * skew(e2);
M156 = -skew(SigB(1:3,5)) + SigB(4:6,1:3) * skew(e2);
M165 = -skew(SigB(4:6,6)) + SigB(4:6,4:6) * skew(e3);
M166 = -skew(SigB(1:3,6)) + SigB(4:6,1:3) * skew(e3);

M = [M;
    zeros(3,12),  M55,  M56;
    zeros(3,12),  M65,  M66;
    zeros(3,12),  M75,  M76;
    zeros(3,12),  M85,  M86;
    zeros(3,12),  M95,  M96;
    zeros(3,12), M105, M106;
    zeros(3,12), M115, M116;
    zeros(3,12), M125, M126;
    zeros(3,12), M135, M136;
    zeros(3,12), M145, M146;
    zeros(3,12), M155, M156;
    zeros(3,12), M165, M166];

RHS2 = SE3Adinv(Z) * SigC * SE3Adinv(Z)' - SigB;
RHS2 = reshape(RHS2, 3, 12);

b = [b; RHS2(:)];
end

function [M, b] = MbMat_2(C,Z,Binv,Yinv,A,X,SigB,SigA,B)
% Construction M and b matrices
e1 = [1;0;0]; e2 = [0;1;0]; e3 = [0;0;1];
Binv = inv(B);
SigBinv = SE3Ad(B) * SigB * SE3Ad(B)';

%% CZB^{-1} = Y^{-1}AX
% Rotation part
M11 = Yinv(1:3,1:3) * A(1:3,1:3) * X(1:3,1:3) * skew(e1);
M13 = -skew(Yinv(1:3,1:3) * A(1:3,1:3) * X(1:3,1:3) * e1);
M15 = -C(1:3,1:3) * Z(1:3,1:3) * skew(Binv(1:3,1:3)*e1);

M21 = Yinv(1:3,1:3) * A(1:3,1:3) * X(1:3,1:3) * skew(e2);
M23 = -skew(Yinv(1:3,1:3) * A(1:3,1:3) * X(1:3,1:3) * e2);
M25 = -C(1:3,1:3) * Z(1:3,1:3) * skew(Binv(1:3,1:3)*e2);

M31 = Yinv(1:3,1:3) * A(1:3,1:3) * X(1:3,1:3) * skew(e3);
M33 = -skew(Yinv(1:3,1:3) * A(1:3,1:3) * X(1:3,1:3) * e3);
M35 = -C(1:3,1:3) * Z(1:3,1:3) * skew(Binv(1:3,1:3)*e3);

% Translation part
M42 = -Yinv(1:3,1:3) * A(1:3,1:3) * X(1:3,1:3);
M43 = -skew(Yinv(1:3,1:3)*A(1:3,1:3)*X(1:3,4) + Yinv(1:3,1:3)*A(1:3,4) + Yinv(1:3,4));
M44 = eye(3);
M45 = -C(1:3,1:3) * Z(1:3,1:3) * skew(Binv(1:3,4));
M46 = C(1:3,1:3) * Z(1:3,1:3);

M = [       M11, zeros(3,3), M13, zeros(3,3), M15, zeros(3,3);
            M21, zeros(3,3), M23, zeros(3,3), M25, zeros(3,3);
            M31, zeros(3,3), M33, zeros(3,3), M35, zeros(3,3);
     zeros(3,3),        M42, M43,        M44, M45,       M46];

% RHS
RHS = - C * Z * Binv + Yinv * A * X;

b = [RHS(1:3,1); RHS(1:3,2); RHS(1:3,3); RHS(1:3,4)];

%% SigBi^{-1} = Ad^{-1}(X) * SigAi * Ad^{-T}(X)
% First block
M51 = -skew(SigBinv(1:3,1)) + SigBinv(1:3,1:3) * skew(e1);
M52 = zeros(3,3);
M61 = -skew(SigBinv(1:3,2)) + SigBinv(1:3,1:3) * skew(e2);
M62 = zeros(3,3);
M71 = -skew(SigBinv(1:3,3)) + SigBinv(1:3,1:3) * skew(e3);
M72 = zeros(3,3);

% Second block
M81 = -skew(SigBinv(1:3,4)) + SigBinv(1:3,4:6) * skew(e1);
M82 = SigBinv(1:3,1:3) * skew(e1);
M91 = -skew(SigBinv(1:3,5)) + SigBinv(1:3,4:6) * skew(e2);
M92 = SigBinv(1:3,1:3) * skew(e2);
M101 = -skew(SigBinv(1:3,6)) + SigBinv(1:3,4:6) * skew(e3);
M102 = SigBinv(1:3,1:3) * skew(e3);

% Third block
M111 = -skew(SigBinv(4:6,1)) + SigBinv(4:6,1:3) * skew(e1);
M112 = -skew(SigBinv(1:3,1));
M121 = -skew(SigBinv(4:6,2)) + SigBinv(4:6,1:3) * skew(e2);
M122 = -skew(SigBinv(1:3,2));
M131 = -skew(SigBinv(4:6,3)) + SigBinv(4:6,1:3) * skew(e3);
M132 = -skew(SigBinv(1:3,3));

% Fourth block
M141 = -skew(SigBinv(4:6,4)) + SigBinv(4:6,4:6) * skew(e1);
M142 = -skew(SigBinv(1:3,4)) + SigBinv(4:6,1:3) * skew(e1);
M151 = -skew(SigBinv(4:6,5)) + SigBinv(4:6,4:6) * skew(e2);
M152 = -skew(SigBinv(1:3,5)) + SigBinv(4:6,1:3) * skew(e2);
M161 = -skew(SigBinv(4:6,6)) + SigBinv(4:6,4:6) * skew(e3);
M162 = -skew(SigBinv(1:3,6)) + SigBinv(4:6,1:3) * skew(e3);

M = [M;
     M51,  M52, zeros(3,12);
     M61,  M62, zeros(3,12);
     M71,  M72, zeros(3,12);
     M81,  M82, zeros(3,12);
     M91,  M92, zeros(3,12);
    M101, M102, zeros(3,12);
    M111, M112, zeros(3,12);
    M121, M122, zeros(3,12);
    M131, M132, zeros(3,12);
    M141, M142, zeros(3,12);
    M151, M152, zeros(3,12);
    M161, M162, zeros(3,12)];

RHS2 = SE3Adinv(X) * SigA * SE3Adinv(X)' - SigBinv;
RHS2 = reshape(RHS2, 3, 12);

b = [b; RHS2(:)];
end

function invX = SE3inv(X)
invX = [X(1:3,1:3)', -X(1:3,1:3)'*X(1:3,4); 0 0 0 1];
end

function A = SE3Adinv(X)
R = X(1:3,1:3);
t = X(1:3,4);

A = [R', zeros(3,3); -skew(R'*t)*R', R'];
end

function A = SE3Ad(X)
R = X(1:3,1:3);
t = X(1:3,4);

A = [R, zeros(3,3); skew(t)*R, R];
end
