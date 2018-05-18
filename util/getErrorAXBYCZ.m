function xyzError = getErrorAXBYCZ(X_f, Y_f, Z_f, XActual, YActual, ZActual)
% Compute the rotational and translational errors for X, Y and Z

xyzError(1) = roterror(X_f, XActual);
xyzError(2) = roterror(Y_f, YActual);
xyzError(3) = roterror(Z_f, ZActual);

xyzError(4) = tranerror(X_f, XActual);
xyzError(5) = tranerror(Y_f, YActual);
xyzError(6) = tranerror(Z_f, ZActual);

end