function [A1, B1, C1, A2, B2, C2, A3, B3, C3] = ...
            generateSetsOfABC(Num, optPDF, Mean, Cov, XActual, YActual, ZActual)
%% Generate sets of {A,B,C} triples 
% Generate three sets of {A,B,C} triples given options for different
% distributions

% Generate constant A1, free B1 and C1
opt = 1;
[A1, B1, C1] =  generateABC(Num, opt, optPDF, Mean, Cov, ...
                            XActual, YActual, ZActual);

% Generate constant C2, free A2 and B2
opt = 3;
[A2, B2, C2] =  generateABC(Num, opt, optPDF, Mean, Cov, ...
                            XActual, YActual, ZActual);

% Generate constant B3, free A3 and C3
opt = 2;
[A3, B3, C3] =  generateABC(Num, opt, optPDF, Mean, Cov, ...
                            XActual, YActual, ZActual);

end

