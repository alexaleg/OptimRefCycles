
function [c,ceq]=f(w)
% This is file flashPTsrk.m
%    Flash at given pressure and temperature
%    w = internal state vector of length 3*NC-1
% Usage for a case with NC=5 components (note: give compositions only for NC-1 components):
%    x0=Initial estimate x, y0= initial estimate y; s0=slack variable based on imaginary liquid state; VF0= initial estimate vf;  w0= [x0, y0, s0 VF0];
%    w=fmincon('1',w0,[],[],[],[],[],[],@flashPTsrk);
%    x =  [w(1:NC-1) 1-sum(w(1:NC-1))]; y =[w(NC:2*NC-2) 1-sum(w(NC:2*NC-2))]; s =[w(2*NC-1:3*NC-2) 1-sum(w(2*NC-1:3*NC-2))]; VF=w(3*NC-1);
%   Vapor_fraction=VF, Liquid_Composition = x, Vapor_Composition = y
% Definition of global variables for recursive routine
global T p zf
% All in SI units
% T = Temperature [K]
% p = Pressure [N/m^2]
% zf = Initial composition

%Number of components
NC=5;

liquid=1;
vapor=2;
% Extract states:
x(1:NC-1)=w(1:NC-1);                   % liquid composition
y(1:NC-1)=w(NC:2*NC-2);                % vapor composition
s(1:NC)=w(2*NC-1:3*NC-2);
VF =w(3*NC-1);                         % vapor fraction
beta(1:NC)=w(3*NC:4*NC-1);
% Composition of last component:
x(NC)=1-sum(x(1:NC-1));
y(NC)=1-sum(y(1:NC-1));


% Fugacity, enthalpy and volumes from SRK
[~,phil,~,~]=srks(x,T,p,liquid);       % liquid
[~,phig,~,~]=srks(y,T,p,vapor);        % vapor
M =1;
K = phil./phig;

f1 = y - beta.*K.*x;                                                % = 0   Algebraic (NC): VLE
f2(1:NC-1) = zf(1:NC-1) - VF.*y(1:NC-1) - (1-VF).*x(1:NC-1);  % = 0   Algebraic(NC-1): Component mass balances

f3(1:NC)=beta-ones-s;
f5(1:NC)=-s;                                                   % < 0   Slack variables
f6 = ones - VF*ones;

c=[f3';f5'; f6']; ceq=[f2'; f1'];
