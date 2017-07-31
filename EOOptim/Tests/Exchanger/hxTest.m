%% Equations oriented HA-02 test
clear
clc

load('../../../ModRun/Modres.mat')

NC = 5;
w0 = [];
for i = [8 11 13]
  w0 = [w0,...
        Result.x(i,1:4), ...
        Result.y(i,1:4), ...
        Result.VF(i), ...
        Result.Zl(i), Result.Zv(i), ...
        Result.sl(i), Result.sv(i), ...
        Result.beta(i), Result.T(i)];
end
Split = 0.65;
QHA_04 = Result.QHA_04*1.0;
h12 = Result.h(12);
w0 = [w0, Split, QHA_04];

lb = zeros(size(w0));
ub = ones(size(w0));

for i = 1:3
  ub(i*3*NC-1) = Inf;
  ub(i*3*NC) = Inf;
end
% ub(end-1)   = Inf;
% lb(end) = -Inf;
ub(end) = Inf;

rho = 2e4;
Plow = 3e5;
flow = Result.n(8);

% Relaxed from 1e-6 to 1e-5 for both optimality and constraint tolerances
options = optimoptions('fmincon','Algorithm','interior-point', 'Display', 'iter-detailed',...
                       'ScaleProblem', 'obj-and-constr', 'OptimalityTolerance', 1e-5,...
                       'ConstraintTolerance', 1e-5);

[w,fval,exitflag,output]=fmincon(@(w)hxObjFunc(w,rho),w0,[],[],[],[],lb,ub,...
                                @(w)hxModel(w, Plow, flow),options);
