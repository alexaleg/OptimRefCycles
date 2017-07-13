function [ Tout_real ,W,hout,hin,Tout] = compressor( Pin,Tin,Pout,x,N)
%This function corresponds to the recursive routine to calculate the
%compressor. This function uses a simplified isoenthropic model.
% % Input
% Pin = Inlet pressure [Pa]
% Tin = Inlet temperature [K]
% Pout = Outlet pressure [Pa]
% x = Inlet molar composition
% N = Inlet molar flow [Kmol/s]
% % Output
% Tout1 = Outlet temperature [K]
% W = Compressor work [Kj/s]
% hout = Outlet molar enthalpy [Kj/Kmol]
% hin = Inlet molar enthalpy [Kj/Kmol]

%% Curves calculation
p = [0.0015   -0.0314    0.1911    0.4217]; %3rd degree polynomial.
rel = Pout/Pin; % Values from 2 to 10
eta = polyval(p,rel);

%% Ideal gas stimation of output temperature
c=cpav(Tin,x); % Calculates the average Cp [Kj/K] for a given mixture at a especified temperature.
R = 8.314; % Universal gas constant [Kj/KmolK]
%eta_0=0.74; % Compressor adiabatic efficiency
k = c/(c-R);
p = (k-1)/k;
%Tout=Tin*(1+((Pout/Pin)^(p)-1)/eta); % Output temperature isoentropic
T0=Tin*((Pout/Pin)^(p/eta)); % Output temperature polytropic

%% Isoentropic calculation

[~,~,Hi,Vi,Si]= srk(x,Tin,Pin,2);                   % Inlet enthalpy and entropy
Tout = fzero(@(T)Sentropy( Si, x ,T, Pout ),T0);     % Output isoentropic temperature
[~,~,Ho,Vo,~]= srk(x,Tout,Pout,2);                  % Output isoentropic enthalpy



%% Real calculations

W=(Ho-Hi)/eta;    % Work calculation for
Ho_real = Hi + W; % Real output enthalpy

Tout_real = fzero(@(T)Centhalpy( Ho_real, x ,T, Pout ),T0); % Real output temperature

W = W * N; % Work for the defined molar flow

% Enthalpies extraction
hout=Ho_real;
hin=Hi;
end
