function [c,ceq]=f(w, z, T, P)

% TESTS
% Passed test for close to vapor

% TODO
% Liquid test
NC = 5;
x = zeros(1, NC);
y = zeros(1, NC);

% Extract states:
x(1:NC-1)=w(1:NC-1);                   % liquid composition
y(1:NC-1)=w(NC:2*NC-2);                % vapor composition
VF =w(2*NC-1);                         % vapor fraction
sl = w(2*NC);                   % liquid slack
sv =w(2*NC+1);                    % vapor slack
beta =w(2*NC+2);                  % relaxation parameter

% Last composition
x(NC)= 1-sum(x);
y(NC)= 1-sum(y);

% Fugacity, and enthalpy from SRK
liquid=1;vapor=2;
[~,phil,hf,~] = srk(x,T,P,liquid);       % liquid
[~,phig,hg,~] = srk(y,T,P,vapor);        % vapor

K = phil./phig;

% Equality constraints
f1 = y - beta.*K.*x;                                          % = 0   Algebraic (NC): VLE
f2 = z - VF.*y - (1-VF).*x;  % = 0   Algebraic(NC-1): Component mass balances
f3 = sum(y) - sum(x);

% Inequality constraints
f4 = 1 - beta - sl; % < 0  liquid slack relaxation
f5 = beta - 1 - sv; % < 0  vapor slack relaxation
% slacks > 0
f6 = - sl;
f7 = - sv;
f8 = -VF*sv;
f9 = -(1 - VF)*sl ;


ceq=[f1'; f2'; f3'; f8'; f9'];
c=[f4'; f5'; f6'; f7'];

% % Include energy balance in case PH flash
% if PH == True
%   fh = hf - VF*hg - (1-VF)*hl;
%   ceq = [ceq; fh'];
