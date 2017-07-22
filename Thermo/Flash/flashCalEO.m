function [x, y, VF, T, Zl, Zv, sl, sv, beta, h] = flashCal(z, P, x0, y0, VF0, T, Zl, Zv, sl, sv, beta, rho, kind, varargin)
%% This is flashCalEO.m - Wrapper function for flash calculations
% Information on the constraints can be found in:
% PTCal.m for PT flash
% PHCal.m for PH flash
% objective.m for objective function

% Input:
%   - z:  molar feed composition (NC)
%   - P: pressure [Pa]
%   - x0: molar composition estimate for liquid phase (NC-1)
%   - y0: molar composition estimate for vapor phase (NC-1)
%   - VF0: estimate vapor fraction
%   - T: temperature(PT flash) / estimate (PH flash)
%   - Zl: liquid phase compressibility estimate
%   - Zv: vapor phase compressibility estimate
%   - sl: liquid phase slack estimate
%   - sv: vapor phase slack estimate
%   - beta: relaxation parameter for VLE
%   - rho: objective coefficient
%   - kind: calculation type 'PT' or 'PH'
%   - varargin: variable input. Not existent for PT. h[kj/kmol] for PH.
% Output:
%   The output is the final iteration value for the same variables.
global flash
flash = flash+1
NC=5; % Number of components


% Algorithm options
options = optimoptions('fmincon','Algorithm','interior-point', 'Display', 'notify-detailed');

if kind == 'PT' % PT flash calculations
  w0= [x0, y0, VF0, Zl, Zv, sl, sv, beta];
  % Boundaries
  lb = zeros(size(w0));
  ub = ones(size(w0));
  lb = zeros(size(w0));
  ub = ones(size(w0));
  ub(2*NC+2:end)=inf;
  [w,fval,exitflag,output]=fmincon(@(w)objective(w,rho),w0,[],[],[],[],lb,ub,@(w)PTCal(w, z, P, T),options);

elseif kind == 'PH' % PT flash calculations
  h = varargin{1}; % Extract h
  w0 = [x0, y0, VF0, Zl, Zv, sl, sv, beta, T];
  % Boundaries
  lb = zeros(size(w0));
  ub = ones(size(w0));
  lb = zeros(size(w0));
  ub = ones(size(w0));
  ub(2*NC+2:end)=inf;
  [w,fval,exitflag,output]=fmincon(@(w)objective(w,rho),w0,[],[],[],[],lb,ub,@(w)PHCal(w, z, P, h),options);
  T = w(end);
else
  disp('Poorly specified flash')
end

%% Data extraction
% Liquid and vapor compositions
x = [w(1:NC-1) 1-sum(w(1:4))];
y = [w(NC:2*NC-2) 1-sum(w(5:8))];

VF = w(2*NC-1);                         % vapor fraction
Zl = w(2*NC);                           % liquid compressibility
Zv = w(2*NC+1);                         % vapor compressibility
sl = w(2*NC+2);                         % liquid slack (imaginary phase)
sv = w(2*NC+3);                         % vapor slack (imaginary phase)
beta = w(2*NC+4);                       % relaxation parameter

% Enthalpy
[~,hf] = srkeo(x,T,P,Zl);       % liquid
[~,hg] = srkeo(y,T,P,Zv);       % vapor
h = VF*hg + (1-VF)*hf;
end
