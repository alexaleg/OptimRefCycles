function [J]=f(w, rho)
%% Objective function for flash calculations
% Assures complementarity of slack phases
%   J = rho*(VF*sv+(1-VF)*sl)
%     if VF = 0, then sl = 0
%     if VF = 1, then sv = 0
%     if 0 < Vf < 1, then sl = sv = 0

% Input:
%   - w: vector of states for optimization
%   - rho: penalty coefficient (supplied by user in wrapper)
% Output:
%   - J: rho*(VF*sv+(1-VF)*sl)

NC = 5;

% Extract states:
sl = w(2*NC+2);                   % liquid slack
sv = w(2*NC+3);                   % vapor slack
VF = w(2*NC-1);

% Penalty function
J = rho*(VF*sv+(1-VF)*sl);
end
