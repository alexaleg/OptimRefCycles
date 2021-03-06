%% Equations oriented compressor test
clear
clc

load('../../../ModRun/Modres.mat')
i = 1;
w0 = [];
w0 = [w0,...
      Result.x(i,1:4), ...
      Result.y(i,1:4), ...
      Result.VF(i), ...
      Result.Zl(i), Result.Zv(i), ...
      Result.sl(i), Result.sv(i), ...
      Result.beta(i), Result.T(i)];

Tout = 360;
w0 = [w0 Tout];

% Relaxed from 1e-6 to 1e-5 for both optimality and constraint tolerances
options = optimoptions('fmincon','Algorithm','interior-point', 'Display', 'iter-detailed',...
                       'ScaleProblem', 'obj-and-constr', 'OptimalityTolerance', 1e-5,...
                       'ConstraintTolerance', 1e-5);

% Manually determined parameter rho
rho = 5e6;
Plow = 3e5;
Phigh = 20e5;
lb = zeros(size(w0));
ub = ones(size(w0));
ub(14:end) = inf;

[w,fval,exitflag,output]=fmincon(@(w)compEOObjFunc(w,rho),w0,[],[],[],[],lb,ub,...
                                 @(w)compEOModel(w, Plow, Phigh),options);
x = [w(1:4) 1-sum(w(1:4))]
y = [w(5:8) 1-sum(w(5:8))]
VF = w(9)
T = w(15)
