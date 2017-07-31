%% Equations oriented stream 8 test
clear
clc

load('../../../ModRun/Modres.mat')

NC = 5;
w0 = [];
i = 8;
w0 = [w0,...
      Result.x(i,1:4), ...
      Result.y(i,1:4), ...
      Result.VF(i), ...
      Result.Zl(i), Result.Zv(i), ...
      Result.sl(i), Result.sv(i), ...
      Result.beta(i), Result.T(i)];

QHA_04 = Result.QHA_04*1.2;
w0 = [w0 QHA_04];

lb = zeros(size(w0));
ub = ones(size(w0));

ub(3*NC-1) = Inf;
ub(3*NC) = Inf;
ub(end) = Inf;
rho = 2e4;
flow = Result.n(8);

% Relaxed from 1e-6 to 1e-5 for both optimality and constraint tolerances
options = optimoptions('fmincon','Algorithm','interior-point', 'Display', 'iter-detailed',...
                       'ScaleProblem', 'obj-and-constr', 'OptimalityTolerance', 1e-5,...
                       'ConstraintTolerance', 1e-5);

[w,fval,exitflag,output]=fmincon(@(w)ObjFunc8(w,rho),w0,[],[],[],[],lb,ub,...
                                @(w)Model8(w, flow),options);
