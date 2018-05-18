function [X, Y, Z] = InitializeXYZ(opt)
% Generate one triple of [X, Y, Z] as the ground truth
if opt == 1
    
    x = randn(6,1); x = x./norm(x); X = expm(se3_vec(x));    % Generate a Random X
    y = randn(6,1); y = y./norm(y); Y = expm(se3_vec(y));    % Generate a Random Y
    z = randn(6,1); z = z./norm(z); Z = expm(se3_vec(z));    % Generate a Random Y

elseif opt == 2
    
    X = [-0.9765   0.0636  -0.2059   0.0215;
         -0.0947  -0.9849   0.1447  -0.0029;
         -0.1936   0.1608   0.9678  -0.0597;
          0        0        0        1.0000;];    
    
    Y = [-0.99908 -0.03266  0.02786  164.226/1000;
          0.02737  0.01553  0.99950  301.638/1000;
         -0.03308  0.99935 -0.01462 -962.841/1000;
          0.00000  0.00000  0.00000  1.00000;];
    
    Z = [0.70063   -0.40451    0.58779  0.006;
         0.69084    0.17849   -0.70063  0.030;
         0.17849    0.89695    0.40451  0.921;
         0.00000    0.00000    0.00000  1.000;];

elseif opt == 3
    
    X(1:3,1:3) = rotx(pi/3)*rotz(pi/4);
%     XActual(1:3,1:3) = rotz(pi/4)*roty(pi/6);
%     YActual(1:3,1:3) = rotz(pi/4)*roty(pi/6)*rotz(pi/5);
    Y(1:3,1:3) = rotz(pi/4)*roty(pi/6);
    Z(1:3,1:3) = rotx(pi/2)*rotz(pi/3)*roty(pi/4);
    
else
    
    fprintf('The given opt = %d is not an option. \n', opt);
    return
end


end
