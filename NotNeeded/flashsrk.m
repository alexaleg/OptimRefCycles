function [cons,f_ine,h,penal] = flashsrk(type,z,x,y,T,P,VF,Zl,Zg,varargin)
% Flash calculations function. This function creates the vector of
% constraints for different possibilities of flash calculation.

% Inputs:
% Type: defines the type of flash calculation TP/PH
% Slack: defines the use of slack variables
% z: total composition of the stream
% x: liquid composition
% y: vapour composition
% T: temperature [K]
% P: pressure [bar]
% VF: vapour fraction
% Zl: liquid compressibility factor
% Zg: gas compressibility factor
% varargin: variable input which can be:
%  for PH flash
%     h: enthalpy [kJ/kmol]
%  for Slacks
%   s: slack (imaginary liquid composition)
%   beta: VLE relaxation factor

% Outputs:
% cons:list of equality constraints
% f_ine: list of inequality constraints(if slacks are on)
% h: enthalpy [kJ/kmol] for TP flash
% penal: Penalty constraint

% Initialize parameters
n_e = 1; % Number of equality constraints
n_i = 1; % Number of inequality constraints
NC = 5;  % Number of components
penal = [];

%% Constraints regardless of the type of flash and slacks
[phil,hl,~,errl,~]=srkseo(x,T,P,Zl);       % liquid
[phig,hg,~,errg,~]=srkseo(y,T,P,Zg);       % vapour

K = phil./phig; %Calculate fugacity

f{n_e}(1:NC-1) = z(1:NC-1) - VF.*y(1:NC-1) - (1-VF).*x(1:NC-1); n_e=n_e+1;   % = 0   Algebraic(NC-1): Component mass balances
f{n_e} = [abs(errl) abs(errg)];                                 n_e=n_e+1;   % = 0   Solution of CEOS

%% Energy balance: Uses it as a constraint for PH flash, claculates enthalpy for TP flash

if type == 'PH' %For PH flash the energy balance becomes a
    a=1;            % Initialize index for variable inputs
    h=varargin{a};  % Extract h from the variable input
    f{n_e} = (h - VF*hg - (1-VF)*hl)*10^-3;                             n_e=n_e+1;   % = 0   Algebraic (1): Energy balance;

elseif type == 'TP' % For TP flash the energy balance is calculated based on both the liquid and gas enthalpies
    a=0;    %Initialize index for variable inputs
    [~,~,hl,~]=srks(x,T,P,1);       % liquid
    [~,~,hg,~]=srks(y,T,P,2);       % vapour
    h=hg*VF+hl*(1-VF); % Value of enthalpy for TP flash

else
    disp('Error specifying flash')
end

%% Slack variables: relaxes the VLE constraints when slack variables are present
% For calculations that involve only one phase or are very close to it, it is necessary to relax the constraints
% The relaxation for this problem is based only on an imaginary liquid
% phase of composition [s]. The VLE calculations are relaxed by a factor
% [beta] which allows the VLE constraints to be continous at all points.
% The inequality constraint allows to carry out this relaxation and by
% meeting the slack constraint the value of beta becomes 1, thus meeting
% the VLE constraint.

if slack == 'on '
    s=varargin{a+1};    % Extract s from the variable input
    beta=varargin{a+2}; % Extract beta from the variable input
    f{n_e} = y - beta.*K.*x;
    n_e=n_e+1;   % = 0   Algebraic Relaxed(NC): VLE

    f_in{n_i}(1:NC)=beta-ones-s;                               n_i=n_i+1;   % < 0   Relaxation
    % Complementary constraints
    f_in{n_i}(1:NC)=-s;                                        n_i=n_i+1;   % = 0   Slack variables
    f_in{n_i}(1:NC)=VF-1;                                      n_i=n_i+1;
    penal = f_in{n_i-2}*(f_in{n_i-1})';


elseif slack == 'off' % If slack variables aren't needed the VLE is directly calculated
    f{n_e} = y - K.*x;                                          n_e=n_e+1;   % = 0   Algebraic (NC): VLE
    f_in=[0];           % Initialize inequality vector

else
    disp('Error specifying slacks')
end

%% Extract vectors of constrains
% The constraints are transposed in order for them to be used by fmincon
% without further modification.
% The order is of them is the same as their order of appearance throughout the code.

cons=[];
n_e;
for i = 1:n_e-1
    cons=[cons;f{i}'];
end

f_ine=[];
for i = 1:n_i-1
    f_ine=[f_ine;f_in{i}'];
end

end
