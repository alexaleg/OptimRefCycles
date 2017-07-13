%% SRK Package
function [phi,H,V,S, errZ, derZ, secDerZ]= srkeo(x,T,P, Z)
% srkeo.m
% Calculations using SRK for a equations oriented solution
% Input:
%     - compositon x (vector)
%     - Temperature T [K]
%     - Pressure P [Pa]
%     - Compressibility factor Z
% Output:
%     - Fugacity coeff.(vector) phi [-]
%     - Enthalpy H [kJ/kmol]
%     - Molar Volume V [m3/mol]
%     - Entropy S [kJ/kmolK]
%     - CEOS errZ
%     - First derivative CEOS derZl
%     - Second derivative CEOS secDerZ

% Based on Linhart & Skogestad, 2008
% Corrected 02, 2016

%% Initialize
%Components N2 C1 C2 C3 C4

NC=5;  % Number of components (1)Nitrogen - (2)Methane - (3)Ethane - (4)Propane - (5)Butane

%Component data from: Reid RC, Prausnitz JM, Poling BE. The properties of gases \& liquids (5th edition). New York: McGraw-Hill, Inc., 2001.

compData.Pc = [33.98  45.99  48.72  42.48  37.96 ]*1.e5;    % [N/m2]
compData.Tc = [126.20 190.56 305.32 369.83 425.12];         % [K]
compData.w  = [0.037  0.011  0.099  0.152  0.200 ];         % acentric factor [-]

R = 8.314; %kJ/kmolK

% Ideal has heat capacity (See Hid below):
compData.Cp(:, 1) = [3.539 -0.261e-3 0.007e-5   0.157e-8 -0.099e-11]' * R;    % Nitrogen
compData.Cp(:, 2) = [4.568 -8.975e-3 3.631e-5  -3.407e-8  1.091e-11]' * R;    % Methane
compData.Cp(:, 3) = [4.178 -4.427e-3 5.660e-5  -6.651e-8  2.487e-11]' * R;    % Ethane
compData.Cp(:, 4) = [3.847  5.131e-3 6.011e-5  -7.893e-8  3.079e-11]' * R;    % Propane
compData.Cp(:, 5) = [5.547  5.536e-3 8.057e-5 -10.571e-8  4.134e-11]' * R;    % n-butane
compData.Tref = 298.15;  % for Ideal gas heat capacity
compData.Pref = 1.e5; % [N/m2]

kinteraction=zeros(NC,NC);      % Here: SRK binary interaction parameters set to zero


x = x';

Pc=compData.Pc;
Tc=compData.Tc;
w=compData.w;
ZRA=0.29056-0.08775*w;
Cp=compData.Cp;
Tref=compData.Tref;
Pref=compData.Pref;

%% Calculations for given T, P  and composition (x)

% All equations nr from Soave,1972
Tre=T./Tc;
Pre=P./Pc;
m=0.480+1.574.*w-0.176.*w.^2; %Eq 15
a=(1+m.*(1-Tre.^0.5)).^2; %Eq 13
Ap=0.42747.*a.*Pre./Tre.^2; %Eq 8
Bp=0.08664.*Pre./Tre; % Eq 9

% Start calculations
% Binary a's:
Ab = (Ap' * Ap).^0.5;

% Mixture a and b
A = x'*(Ab.*(1-kinteraction))*x;
B = Bp*x;

%% Cubic equation of state
% Vectors for polinomial calculations
pol = [1 -1 A-B-B^2 -A*B];
derPol = [3 -2 A-B-B^2];
secDerPol = [6 -2];
% For a given Z calculate the error
errZ = polyval(pol, Z);
derZ = polyval(derPol, Z); % First derivative of the polynomial
secDerZ = polyval(secDerPol, Z); % Second derivative of the polynomial


%% Density (more precisely: molar volume)
% TODO: Vectorize calculation for correction
Vt = Z * R * T / P;
if secDerZ < 0 % Correct liquid SRK-volume using Peneleoux correction
     c=0;
     for i=1:NC
       c=c+x(i) * (0.40768 * (0.29441 - ZRA(i)) * (R * Tc(i)) / (Pc(i))) ;
     end
   V = ((Z * R * T / P)- c);
else  % vapor
   V = Z * R * T / P;
end

%% Fugacity coefficient
% Correction added Nov. 2013:
corrphi = (x'.*Ap.^0.5)*(1-kinteraction);
% factor1 =((corrphi.*2.*Ap.^0.5./A)-Bp./B)
% factor2 = A/B*log((Z+B)/Z)
% Z-B
% term1 = (Z-1).*Bp/B-log(Z-B)
% exponent = (term1-factor2.*factor1)
% exp(exponent)
phi=exp((Z-1).*Bp/B-log(Z-B)...
     -A/B*log((Z+B)/Z).*((corrphi.*2.*Ap.^0.5./A)-Bp./B)); %Eq 21 Soave,1972
% phi=real(phi);
%% Enthalpy Molar
% Calculations from Reid RC, Prausnitz JM, Poling BE. The properties of gases \& liquids (5th edition). New York: McGraw-Hill, Inc., 2001.
% TODO: Vectorize calculation for dadt
dadT=0;
corrb = (R*T/P);
corra = ((R*T)^2/P);

for i = 1:NC
    for j = 1:NC
        dadT = dadT -R / 2 * sqrt(abs(0.42747 / T)) * x(i) * x(j) * (m(j) * sqrt(abs(Ap(i)*(Tc(j)/P*(T^2)/(Pc(j))*(R^2)))) + m(i) * sqrt(abs(Ap(j) * (Tc(i)/P*(T^2)/(Pc(i))*(R^2))))) ;
    end
end

% Excess enthalpy using SRK
Hsrk = - ((A *corra - T*dadT) / (B * corrb )) * log(Vt / (Vt + B * corrb)) + R * T * (1-Z);
% Ideal enthalpy calculations
Hid = (Cp(1,:) * (T - Tref) + 1/2 * Cp(2,:) * (T^2 - Tref^2) + 1/3 * Cp(3,:) * (T^3 - Tref^3) + 1/4 * Cp(4,:) * (T^4 - Tref^4) + 1/5 * Cp(5,:) * (T^5 - Tref^5))*x;
% Total enthalpy calculation
H = Hid - Hsrk;

%% Entropy
% Calculations from Reid RC, Prausnitz JM, Poling BE. The properties of gases \& liquids (5th edition). New York: McGraw-Hill, Inc., 2001.

% Excess entropy using SRK
Ssrk = (dadT/ ( B * corrb )) * log(Vt / ( Vt + B * corrb)) - R*log(Z * (1 - B * corrb/Vt)); %Prausnitz 2001

% Ideal gas entropy
%   Constant pressure contribution
SidP = (Cp(1,:) * log(T/Tref) + Cp(2,:) * (T - Tref) + 1/2 * Cp(3,:) * (T^2 - Tref^2) + 1/3 * Cp(4,:) * (T^3 - Tref^3) + 1/4 * Cp(5,:) * (T^4 - Tref^4))*x;
%   Constant temperature contribution
SidT = R * log(P/Pref);
%   Mixing contribution
SidM = (log(x))'*x;

S = SidP - SidT - SidM - Ssrk;

end
