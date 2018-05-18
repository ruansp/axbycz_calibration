function [w_hat] = skewlog(R)
%SKEWLOG  calculate the log of a rotation matrix
%
%	[W_HAT] = SKEWLOG(R)
%
%
% See also: TWIST, SKEWEXP, TWISTEXP, TWISTLOG.

% $Id: skewlog.m,v 1.1 2009-03-17 16:40:18 bradleyk Exp $
% Copyright (C) 2005, by Brad Kratochvil

% this algorithm from Murray pg. 414

% Commented out for .mex
% global DebugLevel;
%
% if isempty(DebugLevel) || (DebugLevel > 1)
%   if ~isrot(R),
%     error('ROBOTLINKS:skewlog','R is not a rotation matrix')
%   end
% end

if isequalf(R, eye(3), 1e-6)
    % theta = 0;
    w_hat = zeros(3,3);
else
    val = (trace(R)-1)/2;
    % this is to clamp the values to +-1, which sometimes
    % doesn't happen due to floating point noise
    if val > 1
        val = 1;
    elseif val < -1
        val = -1;
    end
    theta = acos(val);
    if 0 == theta
        w_hat = zeros(3,3);
    elseif abs(pi - theta) < 10^-6
        M = (R - eye(3,3))/2;
        m1 = M(1, 1);
        m2 = M(2, 2);
        m3 = M(3, 3);
        w_hat = theta*[               0  -sqrt((m3 - m1 - m2)/2)  sqrt((m2 - m1 - m3)/2);
                 sqrt((m3 - m1 - m2)/2)                      0  -sqrt((m1 - m2 - m3)/2);
                -sqrt((m2 - m1 - m3)/2)  sqrt((m1 - m2 - m3)/2)                      0];
    else
        w_hat = (R-R')/(2*sin(theta))*theta;
    end
end
