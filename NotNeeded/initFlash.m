function [x , y, VF] = initFLash(z, T, P)
% ${2/.*/U/} Description
%  ${1/.*/U/} = ${2/.*/U/}()
%
% Long description

compData.Pc = [33.98  45.99  48.72  42.48  37.96 ]*1.e5;    % [N/m2]
compData.Tc = [126.20 190.56 305.32 369.83 425.12];         % [K]
compData.w  = [0.037  0.011  0.099  0.152  0.200 ];         % acentric factor [-]
Pc= compData.Pc;
Tc = compData.Tc;
w = compData.w;
K = zeros(1,5);
K = exp(5.37.*(1.-w).*(1.-Tc./T)).*(Pc./P);

VF = fzero(@(VF) sum(z./(K+VF)) , 0.5);
if VF < 0
  VF = 0;
elseif VF > 1
  VF = 1;
end

x = z./(1.+VF.*(K-1));
y = K.*x;
end
