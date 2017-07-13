clear all
clc
% % Main Model code
% This code includes a steady state simulation of the whole mini LNG plant.
N_s=17;% Number of streams

%Initialize variables
P  = zeros(1,N_s+1); % Pressures [pa]
Ts = zeros(1,N_s+1); % Stream temperature [K]
n  = ones(1,N_s+1);  % Flows in [kmol/s]
h  = zeros(1,N_s+1); % Molar enthalpy [Kj/kmol]
vf = zeros(1,N_s+1); % Vapor Fraction
H  = zeros(1,N_s+1); % Total Enthalpy [Kj]
x  = zeros(N_s+1,5); % Liquid Composition
y  = zeros(N_s+1,5); % Vapor Composition
z  = zeros(N_s+1,5); % Total Composition
s  = ones(N_s+1,1); % Imaginary liquid slack
beta  = ones(N_s+1,1); % Relaxation parameter

%Set the pressures for the whole system
lp = [7 12 11 13 14 15 16 17 18];
hp = [1 2 3 4 5 6 8 9 10];
P(lp) = 3e5;
P(hp) = 20e5;

%Set known compositions

z(1,:)  = [0.1074 0.3292 0.4096 0.1345 0.0193];
z(2,:)  = [0.1074 0.3292 0.4096 0.1345 0.0193];
z(3,:)  = [0.1074 0.3292 0.4096 0.1345 0.0193];
z(17,:) = [0.1074 0.3292 0.4096 0.1345 0.0193];

%Vapor
y(1,:)  = [0.1074 0.3292 0.4096 0.1345 0.0193];
y(17,:) = [0.1074 0.3292 0.4096 0.1345 0.0193];

% Set known flows
%Load the results from BOG.m
%filename='BOGresults.mat';
%load(filename)
%N = N_ref; %Refrigerant flow [Kmol/s]
N = 0.0434; %Estimate
n = N*n*1;
% Set split
Sf=0.65;

%Set known Temperatures
Ts(17) = -37+273.15;
%Ts(17) = 240;
Ts(2)  = 35+273.15;
Ts(3)  = -25+273.15;
Ts(6)  = 273.15-75;
Ts(8)  = 273.15-146.6;
Ts(10) = Ts(8);
Ts(9)  = Ts(8);
%Ts(14) = 273.15-101.5;

%Set known vapor fractions
vf(17) = 1;
vf(1)  = 1;
vf(2)  = 1;
vf(5)  = 1;
vf(4)  = 0;

%% Compressor
[Ts(1),Wc,h(1),h(17),T1id] = compressor( P(17),Ts(17),P(1),z(17,:),n(1));
H(1)  = n(1)*h(1);
H(17) = n(17)*h(17);

%% Stream calculations
% Vapour only
z(2,:)=z(1,:);
y(2,:)=z(1,:);
[~,~,h(2),~]= srks(z(2,:),Ts(2),P(2),2);
H(2) = h(2)*n(2);
QMR = H(2)-H(1);


%% Flash calculations
QHA_01 = 299.5*1.4;
H(3)=H(2)-QHA_01;
h(3)=H(3)/n(3);
x0=[0.0080 0.0860 0.4590 0.3450];
y0=[0.1620 0.4570 0.3280 0.1840];
VF0=0.65;
T0=Ts(3);
%[ x(3,:), y(3,:),vf(3), Ts(3)] = flashPH(P(3),h(3), x0,y0,z(3,:),T0,VF0);
[ x(3,:), y(3,:),vf(3),~,~] = flashTP(Ts(3),P(3),x0,y0,z(3,:),VF0);
% Assign temperatures
Ts(5)=Ts(3);
Ts(4)=Ts(3);

%Assign flows and compositions
z(5,:)=y(3,:);
y(5,:)=y(3,:);

z(4,:)=x(3,:);
x(4,:)=x(3,:);

n(5)=n(3)*vf(3);
n(4)=n(3)*(1-vf(3));

[~,~,h(5),~] = srks(z(5,:),Ts(5),P(5),2);
[~,~,h(4),~] = srks(z(4,:),Ts(4),P(4),1);

H(5)=n(5)*h(5);
H(4)=n(4)*h(4);


%% Valve FCV-01
%Isoenthalpic Valve
h(7)=h(4);
n(7)=n(4);
z(7,:)=z(4,:);
x0=[0.1 0.1 0.413 0.437]; y0=[0.032 0.296 0.588 0.079]; VF0=0.26; T0=200;
[x(7,:), y(7,:), vf(7), Ts(7)] = flashPH( P(7), h(7), x0, y0, z(7,:), T0, VF0 );
H(7)=n(7)*h(7);


%% Heat exchanger HA-02 first side
%Calculate stream 6
z(6,:)=z(5,:);
n(6)=n(5);
x0=[0.02 0.26 0.59 0.11];y0=[0.28 0.63 0.0 0.001]; VF0=0.5;
[ x(6,:), y(6,:),vf(6), hl, hv] = flashTP(Ts(6),P(6), x0,y0,z(6,:),VF0);
h(6)=hv*vf(6)+hl*(1-vf(6));
H(6)=h(6)*n(6);
QHA_02=H(5)-H(6);

%% Calculate HA-04 first side
%Calculate stream 7
z(8,:)=z(6,:);
n(8)=n(6);
vf(8)=0; % Assumed liquid
x(8,:)=z(6,:);
[Zv,phiv,h(8),Vv]= srks(z(8,:),Ts(8),P(8),1);
H(8)=h(8)*n(8);
QHA_04=H(6)-H(8);

%Split

z(10,:)=z(8,:);
x(10,:)=x(8,:);

z(9,:)=z(8,:);
x(9,:)=x(8,:);

vf(10)=vf(8);
vf(9)=vf(8);

n(10)=n(8)*Sf;
n(9)=n(8)*(1-Sf);

h(10)=h(8);
h(9)=h(8);

H(10)=h(10)*n(10);
H(9)=h(9)*n(9);

Ts(10)=Ts(8);
Ts(9)=Ts(8);

%% FCV-04 // FCV-05
z(11,:)=z(9,:);
n(11)=n(9);
h(11)=h(9);
x0=[0.045 0.5 0.38 0.06]; y0=[0.75 0.24 0.01 0]; VF0=0.1; T0=100;
[x(11,:), y(11,:), vf(11), Ts(11)] = flashPH( P(11), h(11), x0, y0, z(11,:), T0, VF0 );
H(11)=n(11)*h(11);

%Stream 12
z(12,:)=z(10,:);
n(12)=n(10);
h(12)=h(11);
Ts(12)=Ts(11);
x(12,:)=x(11,:);
y(12,:)=y(11,:);
vf(12)=vf(11);
H(12)=n(12)*h(12);

%% HA-04 Second side
n(13)=n(12);
H(13)=H(12)+QHA_04;
h(13)=H(13)/n(13);
z(13,:)=z(12,:);
x0=z(12,1:4); y0=x0; VF0=0.5; T0=200;s0=[ 0 0 0 0 0]; beta=1.5*ones(size(s0));
[x(13,:), y(13,:), vf(13), Ts(13)] = flashPHslack2( P(13), h(13), x0, y0, z(13,:), T0, VF0,s0,beta );

%% QBOG
Tin_bog = -32 + 273.15; %Inlet temperature [K]
m_bog=610; %Gas mass flow [Kg/h] (Neskaa, 2010)
[ QHA_07, n_bog] = BOGC( Tin_bog, m_bog );

%% HA-07 Simplified
%Global energy balance
%This energy is adjusted to the requirements of the natural gas by modifying the refrigerant flow.
%These requirements are in BOG.m.
z(14,:)=z(11,:);
n(14)=n(11);
x0=[0.0027 0.0726 0.7749 0.1435]; y0=[0.2321 0.6364 0.1302 0.0014];
VF0=0.5;
H(14) = H(11) - QHA_07;
h(14) = H(14)/n(14);
T0 = 170;s0=[ 0 0 0 0 0]; beta=1.5*ones(size(s0));
[ x(14,:), y(14,:),vf(14), Ts(14)] = flashPHslack2(P(14),h(14), x0,y0,z(14,:),T0,VF0,s0,beta);

%% %Mixer
n(15)=n(7)+n(13)+n(14);
z(15,:)= (n(7).*z(7,:)+n(13).*z(13,:)+n(14).*z(14,:))/n(15);
H(15)=H(7)+H(13)+H(14);
h(15)=H(15)/n(15);
x0=[0.0014 0.0302 0.5580 0.3548]; y0=[0.1628 0.4856 0.3320 0.0193]; VF0=0.65; T0=195;s0= [ 0 0 0 0 0];beta=1.5*ones(size(s0));
[x(15,:), y(15,:), vf(15), Ts(15)] = flashPH( P(15), h(15), x0, y0, z(15,:), T0, VF0);

%% QHA_02 second side
z(16,:)=z(15,:);
n(16)=n(15);
H(16)=H(15)+QHA_02;
h(16)=H(16)/n(16);
x0=[0.0008 0.0131 0.2491 0.5045]; y0=[0.1134 0.3470 0.4187 0.1136]; VF0=0.90; T0=230;s0= [ 1 1 1 1 1];beta=0.9*ones(size(s0));
[x(16,:), y(16,:), vf(16), Ts(16),s(16,:),beta(16,:)] = flashPHslack2( P(16), h(16), x0, y0, z(16,:), T0, VF0,s0,beta );
%% Q_env
%It is calculated by carrying out an energy balance
Q_env=H(17)-H(16);%This corresponds to the heat from the environment
z(18,:)=z(17,:);
n(18)=n(17);
H(18)=Q_env+H(16);
h(18)=H(18)/n(18);
h(18)= h(17);
x0=[0.0004 0.0038 0.1199 0.3908];
%y0=[0.1088 0.1909 0.4453 0.2144];
%x0 = z(17, 1:4);
y0 = z(17, 1:4);
VF0=1; T0=240; s0= [ 1 1 1 1 1];beta=0.9*ones(size(s0));
[x(18,:), y(18,:), vf(18), Ts(18),s(17,:),beta(17,:)] = flashPHslack2( P(18), h(17), x0, y0, z(18,:), T0, VF0,s0,beta );
%[ x(18,:), y(18,:),vf(18), hl, hv] = flashTPslack(Ts(17),P(18), x0,y0,z(18,:),VF0,s0,beta);
Ts(18)=Ts(17);
%h(18)= (1-vf(18))*hl+vf(18)*hv;
H(18)=n(18)*h(18);
x(18,:)*(1-vf(18))+y(18,:)*vf(18);
%Q_env= H(18)-H(16);

%% Temperature Crosses
Cross=[];
%BOG Temperatures from BOG.m
TBOGin = -32 + 273.15; %Inlet temperature [K]
TBOGout = -154 +273.15; %Outlet temperature [K]
% % HA - 02
Cross(1) = -Ts(16) + Ts(5);
Cross(2) = -Ts(15) + Ts(6);
% % HA - 04
Cross(3) = - Ts(13) + Ts(6);
Cross(4) =  -Ts(12) + Ts(8);

% % HA - 07
Cross(5) = - Ts(14) + TBOGin;
Cross(6) = - Ts(11) + TBOGout;

%% Output information
% Result = struct('z',[],'x',[],'y',[],'P',[],'T',[],'VF',[],'h',[],'n',[]);
%
T=Ts;
VF=vf;
% Result.x=x;
% Result.y=y;
% Result.VF=vf;
% Result.P=P;
% Result.T=Ts
% Result.h=h;
% Result.z=z;
% Result.n=n;
% % Save information to file
filename='SSresults.mat';
save(filename)
