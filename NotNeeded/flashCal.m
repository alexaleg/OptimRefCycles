function [x, y, VF, T] = flashCal(z, P, x0, y0, VF0, T, kind, varargin)

% TODO Update comments

NC=5;

sl = 0;
sv = sl;
beta = 1;

options = optimoptions('fmincon','Algorithm','interior-point','TolCon', 1e-6,'TolFun', 1e-4, 'MaxIter', 10000, 'TolX', 1e-13);
if kind == 'TP'
  w0= [x0, y0, VF0, sl, sv, beta];
  lb=zeros(size(w0));
  ub=ones(size(w0));
  lb=zeros(size(w0));
  ub=ones(size(w0));
  lb(2*NC:end)=-inf;
  ub(2*NC:end)=inf;
  [w,fval,exitflag,output]=fmincon('1',w0,[],[],[],[],lb,ub,@(w)TPCaleo(w, z, T, P),options);
elseif kind == 'PH'
  h = varargin{1};
  w0= [x0, y0, VF0, sl, sv, beta, T];
  lb=zeros(size(w0));
  ub=ones(size(w0));
  lb=zeros(size(w0));
  ub=ones(size(w0));
  lb(2*NC:end-1)=-inf;
  ub(2*NC:end)=inf;
  [w,fval,exitflag,output]=fmincon('1',w0,[],[],[],[],lb,ub,@(w)PHCal(w, z, P, h),options);
  T = w(end);
else
  print('Poorly specified flash')
end
%% Data extraction
%Liquid and vapor compositions
x = [w(1:NC-1) 1-sum(w(1:4))];
y = [w(NC:2*NC-2) 1-sum(w(5:8))];
%Vapor fraction
VF = w(2*NC-1);

end
