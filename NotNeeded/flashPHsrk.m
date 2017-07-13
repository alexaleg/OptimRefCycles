function [c,ceq]=f(w, hf, p, zf)
%    Flash at given pressure and enthalpy
%    w = internal state vector of length 4*NC
% Usage for a case with NC=5 components (note: give composition only for NC-1 component):
%    x0= Initial estimates x, y0= Initial estimates yo; VF0= Initial estimate Vf; T0= Initial estimate T;
%    Temperature=T, Vapor_fraction=VF, Liquid_Composition = x, Vapor_Composition = y
% Global variables for recursive routine


% % Data (parameters and independent variables)
%Number of components
NC=5;
% zf % Specified feed composition
% p % Specified pressure [N/m2]
% hf % Specified enthalpy [kJ/Kmol] (ref: h=0 for ideal gas at 298K)

liquid=1;vapor=2;
R=8.13;                                % kJ/Kmol K

% Extract states:
x(1:NC-1)=w(1:NC-1);                 % Liquid composition
y(1:NC-1)=w(NC:2*NC-2);              % Vapor composition
s(1:NC)=w(2*NC-1:3*NC-2);            % Slack variable liquid
VF =w(3*NC-1);                       % Vapor fraction
T =w(3*NC);                          % Temperature [K]
beta=w(3*NC+1:4*NC);                 % Relaxation parameter

% Composition of last component:
x(NC)=1-sum(x(1:NC-1));
y(NC)=1-sum(y(1:NC-1));

% Fugacity, enthalpy and volumes from SRK
[Zl,phil,hl,Vlm]=srks(x,T,p,liquid);       % liquid
[Zg,phig,hg,Vgm]=srks(y,T,p,vapor);        % vapor

K = phil./phig;

M = 1;
f1 = y - beta.*K.*x;                                          % = 0   Algebraic (NC): VLE
f2(1:NC-1) = zf(1:NC-1) - VF.*y(1:NC-1) - (1-VF).*x(1:NC-1);  % = 0   Algebraic(NC-1): Component mass balances
f3 = hf - VF*hg - (1-VF)*hl;                                  % = 0   Algebraic (1): Energy balance
f4(1:NC)= ones - beta - s;                                    % < 0   Relaxation
f5(1:NC)=-s;                                                  % = 0   Slack variables
f6 = -ones+VF;
%f7 = M*f5*f6';
f7 = M*(f5.*f6);
c=[f4';f5';f6';f7']; ceq=[f2';f1';f3'];
