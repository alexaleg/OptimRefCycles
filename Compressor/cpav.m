function [ cpa ] = cpav( T,x )
% % Calculates the average Cp [Kj/K] for a mixture. 
% Given T = Temperature [K]
% x = Mixture molar composition

% % Uses the parameters from srks.m
% Ideal has heat capacity (See Hid below): 
Cp(:, 1) = [7.440 -0.324e-2 6.400e-6 -2.790e-9]'; %Nitrogen
Cp(:, 2) = [4.598 1.245e-2 2.86e-6 -2.703e-9]'; % Methane
Cp(:, 3) = [1.292 4.254e-2 -1.657e-5 2.081e-9]'; % Ethane
Cp(:, 4) = [-1.009 7.315e-2 -3.789e-5 7.678e-9]';    % Propane
Cp(:, 5) = [2.266 7.913e-2 -2.647e-5 -0.674e-9]';    % n-butane
NC=5;
cpa=0;

% Ponderated sum of each components Cp
for i=1:NC 
    cpa = cpa +(x(i)* (Cp(1,i) + Cp(2,i) * T + Cp(3,i) * T^2 +  Cp(4,i) * T^3));
end
cpa= cpa*4.184; % Unit conversion from [Kcal] to [kJ]
end

