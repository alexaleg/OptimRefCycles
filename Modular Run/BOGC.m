function [ QBOG,n_bog ] = BOGC( Tin_bog, m_bog )
% Calculation of the energy requirements of the BOG for a variable inlet
% temperature and flow.

z_bog=[0.11 0.89 0 0 0]; %Gas Mole compositions
M_N= 28; %Nitrogen (N2) molecular mass [Kg/Kmol]
M_C1=16; %Methane (CH4) molecular mass [Kg/Kmol]
n_bog= m_bog/(z_bog(1)*M_N+z_bog(2)*M_C1)*(1/3600); %BOG Molar Flow [Kmol/s]


Tout_bog = -154 +273.15; %Outlet temperature [K]
P_bog = 18.2e5;% Pressure [Pa] (Assumed constant through the exchanger)

[~,~,Hin,~]= srk(z_bog,Tin_bog,P_bog,2); %Inlet SRK Calculations. Initial phase vapor only.
[~,~,Hout,~]= srk(z_bog,Tout_bog,P_bog,1); %Outlet SRK Calculations. Outlet liquid only. 

QBOG= n_bog*(Hout-Hin); %Heat flow from the gas [Kj/s]
end
