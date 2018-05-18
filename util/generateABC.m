function [A, B, C] =  generateABC(length, optFix, optPDF, M, Sig, X, Y, Z)
%% Data generation for AXB = YCZ problem
% Input:
%       length: number of generated data pairs
%       optFix: option for fixing different data streams
%       optPDF: option for generating data using different distributions
%       M:      mean of perturbance in lie algebra
%       Sig:    covariance of perturbance in lie algebra
%       X, Y, Z: ground truths
% Output:
%       A, B, C: 4 x 4 x length or 4 x 4 
%                noise-free data streams with correspondence
%
% Authors: Qianli Ma, qianli.ma622@gmail.com; 
%          Zachariah Goh, zach_goh@yahoo.com

%% Times of simulation steps
len = length; 
% e = precision;
% digits(64);

% a = randn(6,1); a = a./norm(a); A_initial = expm(se3_vec(a));    % Generate a Random A
% b = randn(6,1); b = b./norm(b); B_initial = expm(se3_vec(b));    % Generate a Random A

%using puma560 to generate tranformation A_inittial, B_initial or C_initial
qz1 = [pi/6, pi/3, pi/4,  pi/4, -pi/4, 0];
% qz1 = [pi/3, pi/4, pi/3, -pi/4,  pi/4, 0];
qz2 = [pi/3, pi/4, pi/3, -pi/4,  pi/4, 0];
% qz3 = [pi/3, pi/4, pi/3, -pi/4,  pi/4, 0];
qz3 = [pi/4, pi/3, pi/3,  pi/6, -pi/4, 0];

dataGenMode = 3;
if dataGenMode == 1
    
    % The following instantiation of puma object is time consuming
    mdl_puma560; % include in the puma560 parameters from "rvctools" toolbox
    A_initial = p560.fkine(qz1);
    B_initial = p560.fkine(qz2);
    C_initial = p560.fkine(qz3);
    
elseif dataGenMode == 2
    
    A_initial = ...    
       [ 0.2294   -0.1951   -0.9536   -0.1038;
         0.7098    0.7039    0.0268   -0.2332;
         0.6660   -0.6830    0.3000    0.2818;
              0         0         0    1.0000]; %  p560.fkine(qz1)
    B_initial = ...
       [ 0.0268   -0.7039   -0.7098    0.0714;
        -0.9536    0.1951   -0.2294   -0.1764;
         0.3000    0.6830   -0.6660    0.2132;
              0         0         0    1.0000]; %  p560.fkine(qz2);

    C_initial = ...
       [-0.0335   -0.4356   -0.8995   -0.0128;
         0.4665    0.7891   -0.3995   -0.2250;
         0.8839   -0.4330    0.1768    0.1756;
             0         0         0    1.0000]; %  p560.fkine(qz3);
elseif dataGenMode == 3
    
    a = randn(6,1); a = a./norm(a); A_initial = expm(se3_vec(a));
    b = randn(6,1); b = b./norm(b); B_initial = expm(se3_vec(b));
    c = randn(6,1); c = c./norm(c); C_initial = expm(se3_vec(c));

end

if optFix == 1 % Fix A, randomize B and C
    % This can be applied to both serial-parallel and dual-robot arm
    % calibrations
    
    B = zeros(4, 4, len);
    C = zeros(4, 4, len);
    
    for m = 1:1:len
        
        if optPDF == 1
            B(:,:,m) = expm(se3_vec(mvg(M, Sig, 1)))*B_initial;
        elseif optPDF == 2
            B(:,:,m) = B_initial*expm(se3_vec(mvg(M, Sig, 1)));
        elseif optPDF == 3
            gmean = [0; 0; 0; 0; 0; 0];
            % Assume Sig is a matrix with the same diagonal values
            B(:,:,m) = sensorNoise(B_initial, gmean, Sig(1), 1);
        end
        
        C(:,:,m) = Y \ (A_initial * X * B(:,:,m) / Z);
        A(:,:,m) = A_initial;
    end
    
%     A = A_initial;
    
elseif optFix == 2 % Fix B, randomize A and C
    % This can be applied to both serial-parallel and dual-robot arm
    % calibrations
    A = zeros(4, 4, len);
    C = zeros(4, 4, len);
    
    for m = 1:1:len
        
        if optPDF == 1
            A(:,:,m) = expm(se3_vec(mvg(M, Sig, 1)))*A_initial;
        elseif optPDF == 2
            A(:,:,m) = A_initial*expm(se3_vec(mvg(M, Sig, 1)));
        elseif optPDF == 3
            gmean = [0; 0; 0; 0; 0; 0];
            A(:,:,m) = sensorNoise(A_initial, gmean, Sig(1), 1);
        end
            
        C(:,:,m) = Y \ (A(:,:,m) * X * B_initial / Z);
    end
    
    B = B_initial;
    
elseif optFix == 3 % Fix C, randomize A and B
    % This is only physically achievable on multi-robot hand-eye
    % calibration
    A = zeros(4, 4, len);
    B = zeros(4, 4, len);
    % B_inv = zeros(4, 4, len);

    for m = 1:1:len
        if optPDF == 1
            B(:,:,m) = expm(se3_vec(mvg(M, Sig, 1)))*B_initial;
        elseif optPDF == 2
            B(:,:,m) = B_initial*expm(se3_vec(mvg(M, Sig, 1)));
        elseif optPDF == 3
            gmean = [0; 0; 0; 0; 0; 0];
            B(:,:,m) = sensorNoise(B_initial, gmean, Sig(1), 1);
        end
        %B_inv(:,:,m) = expm(se3_vec(mvg(M, Sig, 1)))*B_initial;
        %B(:,:,m) = inv(B_inv(:,:,m));
        A(:,:,m) = ( Y * C_initial * Z / B(:,:,m) ) / X;
        C(:,:,m) = C_initial;
    end
    
%     C = C_initial;
    
elseif optFix == 4
    % This is for testing tranditional AXBYCZ solver that demands the
    % correspondence between the data pairs {A_i, B_i, C_i}
    A = zeros(4, 4, len);
    B = zeros(4, 4, len);
    C = zeros(4, 4, len);

    for m = 1:1:len
        A(:,:,m) = expm(se3_vec(mvg(M, Sig, 1)))*A_initial;
        C(:,:,m) = expm(se3_vec(mvg(M, Sig, 1)))*C_initial;
        B(:,:,m) = X \ (A(:,:,m) \ Y * C(:,:,m) * Z);
    end    
    
end

end