function [ x, y,VF, hl, hv  ] = flashTPslack(t,P, x0,y0,z,VF0,s0,beta)
% % Recursive routine for TPFlash calculation. 
% Extracts liquid composition, vapor composition, vapour
% fraction, molar enthalpy for liquid and for vapor phase, given Temperature and Pressure. 
%Definition of global variables for recoursive routine 
global T p zf
T=t; p=P; zf=z; 
% zf % Specified feed composition 
% p % Specified pressure [N/m2]
% T % Specified temperature [K]

NC=5;
%Inital estimates for solver
w0= [x0, y0, s0, VF0,beta];    

%s0 initial slack value
%x0 initial estimate for liquid composition
%y0 initial estimate for vapor composition
%vf0 initial estimate for vaporr fraction

%Boundary constrains. All the states must lie between 0 and 1. 
lb=zeros(size(w0));
ub=ones(size(w0));
ub(2*NC-1:3*NC-2)=inf;
lb(3*NC-1)=1;


options = optimoptions('fmincon','Algorithm','sqp','TolCon', 1e-10,'TolFun', 1e-4, 'MaxIter', 3000, 'TolX', 1e-13, 'ScaleProblem','obj-and-constr');
%Solution subject to the equilibrium and solving the EOS. For further
%details check flashTPsrk.m
w=fmincon('1',w0,[],[],[],[],lb,ub,@flashTPsrkslack);

% %Data extraction
%NC=(length(w)+1)/2; 

%Compositions for each phase (Final component calculated)
x =[w(1:NC-1) 1-sum(w(1:NC-1))];
y =[w(NC:2*NC-2) 1-sum(w(NC:2*NC-2))];
VF = w(3*NC-1);

%Compresibility, fugacity, enthalpy and volume calculations
%For the purpose of this project, only enthalpy is needed. 
[Zv,phiv,hv,Vv]= srks(y,T,p,2);
[Zl,phil,hl,Vl]= srks(x,T,p,1);

end

