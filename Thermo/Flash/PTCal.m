function [c,ceq]=f(w, z, P, T)
%% PT flash constraints calculations:
% Input:
%   - w: vector of states for optimization
%   - z: feed composition
%   - P: pressure [Pa]
%   - T: temperature [K]
% Output:
%   - c: nonequality constraints
%      * K value relaxation
%      * slack phases
%      * first derivative condition CEOS
%      * second derivative condicion CEOS
%   - ceq: equality constraints
%      * VLE
%      * component mass balances
%      * sum(mole fractions) = 1
%      * CEOS

%% Initialization
NC = 5;

x = zeros(1, NC);
y = zeros(1, NC);

% Extract states:
x(1:NC-1)=w(1:NC-1);                   % liquid composition
y(1:NC-1)=w(NC:2*NC-2);                % vapor composition
VF =w(2*NC-1);                         % vapor fraction
Zl = w(2*NC);                          % liquid compressibility
Zv = w(2*NC+1);                        % vapor compressibility
sl = w(2*NC+2);                        % liquid slack (imaginary phase)
sv =w(2*NC+3);                         % vapor slack (imaginary phase)
beta =w(2*NC+4);                       % relaxation parameter

% Last composition
x(NC)= 1-sum(x);
y(NC)= 1-sum(y);

%% Calculations
% Fugacity, and enthalpy from SRK

[phil,hf,~,~,fZl, derZl, secDerZl] = srkeo(x,T,P,Zl);       % liquid
[phig,hg,~,~,fZv, derZv, secDerZv] = srkeo(y,T,P,Zv);        % vapor

K = phil./phig;

% Equality constraints
f1 = y - beta.*K.*x;         % = 0   Algebraic (NC): VLE
f2 = z - VF.*y - (1-VF).*x;  % = 0   Algebraic(NC-1): Component mass balances
f3 = sum(y) - sum(x);        % = 0   sum(y) = sum(x) = 1

% Inequality constraints
% Relaxation
f4 = 1 - beta - sl; % < 0  liquid slack relaxation
f5 = beta - 1 - sv; % < 0  vapor slack relaxation

% Compressibility factor
% fZl, fZv = 0; Zl, Zv fulfill the CEOS
M = 5e6; % Relaxation parameter for 2nd derivative
% 1st derivative equal for both phases
f6 = - derZl;
f7 = - derZv;
%  Second derivative
f8 = secDerZl- M*sl;       % f''(Zl) <= M*sl
f9 = - secDerZv - M*sv;    % -f''(Zv) <= M*sv

ceq=[f1'; f2'; f3; fZl; fZv];
c=[f4; f5; f6; f7; f8; f9];
