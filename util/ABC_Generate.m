function [A, B, C] =  ABC_Generate(length, opt, M, Sig, X, Y, Z)
%% Data generation for AXB = YCZ problem

%% Times of simulation steps
len = length; % number of generated data pairs
% e = precision;
% digits(64);

% a = randn(6,1); a = a./norm(a); A_initial = expm(se3_vec(a));    % Generate a Random A
% b = randn(6,1); b = b./norm(b); B_initial = expm(se3_vec(b));    % Generate a Random A

%using puma560 to generate tranformation A_inittial, B_initial or C_initial
qz1 = [pi/6, pi/3, pi/4,  pi/4, -pi/4, 0];
qz2 = [pi/3, pi/4, pi/3, -pi/4,  pi/4, 0];
qz3 = [pi/4, pi/3,-pi/3, -pi/6,  pi/4, 0];

a = randn(6,1); a = a./norm(a); A_initial = expm(se3_vec(a));
b = randn(6,1); b = b./norm(b); B_initial = expm(se3_vec(b));
c = randn(6,1); c = c./norm(c); C_initial = expm(se3_vec(c));


if opt == 1 % Fix A, randomize B and C
    % This can be applied to both serial-parallel and dual-robot arm
    % calibrations
    
    B = zeros(4, 4, len);
    C = zeros(4, 4, len);
    
    for m = 1:1:len
        B(:,:,m) = expm(se3_vec(mvg(M, Sig, 1)))*B_initial;
        C(:,:,m) = Y \ A_initial * X * B(:,:,m) / Z;
    end
    
    A = A_initial;
    
elseif opt == 2 % Fix B, randomize A and C
    % This can be applied to both serial-parallel and dual-robot arm
    % calibrations
    A = zeros(4, 4, len);
    C = zeros(4, 4, len);
    
    for m = 1:1:len
        A(:,:,m) = expm(se3_vec(mvg(M, Sig, 1)))*A_initial;
        C(:,:,m) = Y \ A(:,:,m) * X * B_initial / Z;
    end
    
    B = B_initial;
    
elseif opt == 3 % Fix C, randomize A and B
    % This is only physically achievable on multi-robot hand-eye
    % calibration
    A = zeros(4, 4, len);
    B = zeros(4, 4, len);
    
    for m = 1:1:len
        A(:,:,m) = expm(se3_vec(mvg(M, Sig, 1)))*A_initial;
        % zg: need to have brackets when using mldivide notation
        B(:,:,m) = X \ (A(:,:,m) \ Y * C_initial * Z);
    end
    
    C = C_initial;
    
elseif opt == 4
    % This is for testing tranditional AXBYCZ solver that demands the
    % correspondence between the data pairs {A_i, B_i, C_i}
    A = zeros(4, 4, len);
    B = zeros(4, 4, len);
    C = zeros(4, 4, len);
    
    for m = 1:1:len
        A(:,:,m) = expm(se3_vec(mvg(M, Sig, 1)))*A_initial;
        C(:,:,m) = expm(se3_vec(mvg(M, Sig, 1)))*C_initial;
        % zg: need to have brackets when using mldivide notation
        B(:,:,m) = X \ (A(:,:,m) \ (Y * C(:,:,m) * Z));
        % B(:,:,m) = inv(X) * inv(A(:,:,m)) * Y * C(:,:,m) * Z;
    end
end

end