%% Main Model code
% This code includes a steady state simulation of the whole mini LNG plant.

clear all
clc
% TODO
% Include general comments about how it was initialized
% Clean and organize this code.
% Remove old comments
% Set up results for plantwide optimization

global flash
flash = 0;
N_s = 17;  % Number of streams
rho = 5e4; % Penalty factor

%Initialize variables
P     = zeros(1,N_s+1); % Pressures [pa]
T     = zeros(1,N_s+1); % Stream temperature [K]
n     = ones(1,N_s+1);  % Flows in [kmol/s]
h     = zeros(1,N_s+1); % Molar enthalpy [Kj/kmol]
vf    = zeros(1,N_s+1); % Vapor Fraction
H     = zeros(1,N_s+1); % Total Enthalpy [Kj]
x     = zeros(N_s+1,5); % Liquid Composition [molar]
y     = zeros(N_s+1,5); % Vapor Composition [molar]
z     = zeros(N_s+1,5); % Total Composition [molar]
Zl    = ones(N_s+1,1);  % Compressibility factor liquid
Zv    = ones(N_s+1,1);  % Compressibility factor vapor
sl    = ones(N_s+1,1);  % Imaginary liquid slack
sv    = ones(N_s+1,1);  % Imaginary vapor slack
beta  = ones(N_s+1,1);  % Relaxation parameter

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
N = 0.0434; %Estimate refrigerant flow [Kmol/s]
n = N*n*1;

% Set split
Sf = 0.65;

%Set known Temperatures
T(17) = -37+273.15;
T(2)  = 35+273.15;
T(3)  = -25+273.15;
T(6)  = 273.15-75;
T(8)  = 273.15-146.6;
T(10) = T(8);
T(9)  = T(8);


%Set known vapor fractions
vf(17) = 1;
vf(1)  = 1;
vf(2)  = 1;
vf(5)  = 1;
vf(4)  = 0;

%% Compressor
[T(1),Wc,h(1),h(17), ~] = compressor( P(17),T(17),P(1),z(17,:),n(1));
H(1)  = n(1)*h(1);
H(17) = n(17)*h(17);

% Initial guesses
Zl(1) = 0.2;
Zv(1) = 0.8;
sl(1) = 1;
sv(1) = 0;
beta(1) = 1;

[x(1,:), y(1,:), vf(1), T(1), Zl(1), Zv(1), sl(1), sv(1), beta(1), h(1)] = ...
flashCalEO(z(1,:), P(1), z(1,1:4), z(1,1:4), 1, T(1), Zl(1), Zv(1),...
           sl(1), sv(1), beta(1) , rho, 'PH', h(1));

%% Stream calculations
% Vapour only
n(2) = n(1);
z(2,:) = z(1,:);
y(2,:) = z(1,:);
[Zv(2),~,h(2),~]= srk(z(2,:),T(2),P(2),2);

QMR = -133.1815;
H(2) = H(1) + QMR;
h(2) = H(2) / n(2);
[x(2,:), y(2,:), vf(2), T(2), Zl(2), Zv(2), sl(2), sv(2), beta(2), h(2)] =...
flashCalEO(z(2,:), P(2), x(1,1:4), z(2,1:4), vf(1), T(2), Zl(1), Zv(1) ,...
           sl(1), sv(1), beta(1) , rho, 'PH', h(2));


%% Flash calculations
QHA_01 = 299.5*1.4;
H(3)   = H(2)-QHA_01;
h(3)   = H(3)/n(3);
[x(3,:), y(3,:), vf(3), T(3), Zl(3), Zv(3), sl(3), sv(3), beta(3), h(3)] = ...
flashCalEO(z(3,:), P(3), z(1,1:4), z(1,1:4), 0.8, T(3),...
           Zl(1), Zv(1), sl(1), sv(1), beta(1) , rho, 'PH', h(3));

% Assign temperatures
T(5)=T(3);
T(4)=T(3);

%Assign flows and compositions
z(5,:)=y(3,:);
y(5,:)=y(3,:);

z(4,:)=x(3,:);
x(4,:)=x(3,:);

n(5)=n(3)*vf(3);
n(4)=n(3)*(1-vf(3));

[Zv(5),~,h(5),~] = srk(z(5,:),T(5),P(5),2);
[Zl(4),~,h(4),~] = srk(z(4,:),T(4),P(4),1);

H(5)=n(5)*h(5);
H(4)=n(4)*h(4);

%% Initialization for stream 6
[x(5,:), y(5,:), vf(5), T(5), Zl(5), Zv(5), sl(5), sv(5), beta(5), h(5)] = ...
flashCalEO(z(5,:), P(5), z(5,1:4), z(5,1:4), vf(5), T(5),...
           Zl(3), Zv(5), sl(3), sv(3), beta(3) , rho, 'PH', h(5));

%% Initialization for stream 7
[x(4,:), y(4,:), vf(4), T(4), Zl(4), Zv(4), sl(4), sv(4), beta(4), h(4)] = ...
flashCalEO(z(4,:), P(4), z(4,1:4), z(4,1:4), vf(4), T(4),...
          Zl(4), Zv(3), sl(3), sv(3), beta(3) , rho, 'PH', h(4));


% %% Valve FCV-01
%Isoenthalpic Valve stream 7
h(7)=h(4);
n(7)=n(4);
z(7,:)=z(4,:);

[x(7,:), y(7,:), vf(7), T(7), Zl(7), Zv(7), sl(7), sv(7), beta(7), h(7)] =...
flashCalEO(z(7,:), P(7), x(4,1:4), y(4,1:4), 0.6, T(4), Zl(4), Zv(4),...
           sl(4), sv(4), beta(4) , rho, 'PH', h(7));

H(7)=n(7)*h(7);

%% Heat exchanger HA-02 first side
%Calculate stream 6
z(6,:)=z(5,:);
n(6)=n(5);

[x(6,:), y(6,:), vf(6), T(6), Zl(6), Zv(6), sl(6), sv(6), beta(6), h(6)] =...
flashCalEO(z(6,:), P(6), x(5,1:4), y(5,1:4), 0.5, T(6), Zl(5), Zv(5) ,...
           sl(5), sv(5), beta(5) , rho, 'PT');

H(6)=h(6)*n(6);
QHA_02=H(5)-H(6);

%% Calculate HA-04 first side

%Calculate stream 8
z(8,:) = z(6,:);
n(8)   = n(6);


% Liquid initialization
[Zl(8),~, ~,~] = srk(z(8,:),T(8),P(8),1);
vf(8)  = 0.5;
x(8,:) = z(6,:);

[x(8,:), y(8,:), vf(8), T(8), Zl(8), Zv(8), sl(8), sv(8), beta(8), h(8)] =...
flashCalEO(z(8,:), P(8), z(8,1:4), z(6,1:4), vf(8), T(8), Zl(8), Zv(6),...
           sl(6), sv(6), beta(6) , rho, 'PT');

H(8) = h(8)*n(8);
QHA_04 = H(6)-H(8);

%% Split
z(10,:) = z(8,:);
x(10,:) = x(8,:);

z(9,:) = z(8,:);
x(9,:) = x(8,:);

vf(10) = vf(8);
vf(9)  = vf(8);

n(10) = n(8)*Sf;
n(9)  = n(8)*(1-Sf);

h(10) = h(8);
h(9)  = h(8);

H(10) = h(10)*n(10);
H(9)  = h(9)*n(9);

T(10) = T(8);
T(9)  = T(8);

%% FCV-04 // FCV-05
z(11,:) = z(9,:);
n(11)   = n(9);
h(11)   = h(9);

% Expected values
vf(11)   = 0.3;
sl(11)   = 0;
sv(11)   = 0;
beta(11) = 1;

[x(11,:), y(11,:), vf(11), T(11), Zl(11), Zv(11),...
 sl(11), sv(11), beta(11), h(11)] =...
 flashCalEO(z(11,:), P(11), z(11,1:4), z(11,1:4), vf(11), T(8), Zl(8), Zv(8),...
            sl(11), sv(11), beta(11) , rho, 'PH', h(11));

H(11)=n(11)*h(11);

%Stream 12
z(12,:) = z(10,:);
n(12)   = n(10);
h(12)   = h(11);
T(12)   = T(11);
x(12,:) = x(11,:);
y(12,:) = y(11,:);
vf(12)  = vf(11);
H(12)   = n(12)*h(12);

%% HA-04 Second side
n(13)=n(12);
H(13)=H(12)+QHA_04;
h(13)=H(13)/n(13);
z(13,:)=z(12,:);
% Expected value
vf(13) = 0.5;

[x(13,:), y(13,:), vf(13), T(13), Zl(13), Zv(13),...
 sl(13), sv(13), beta(13), h(13)] =...
 flashCalEO(z(13,:), P(13), x(11,1:4), y(11,1:4), vf(13), T(11) , Zl(11), Zv(11),...
            sl(11), sv(11), beta(11) , rho, 'PH', h(13));



%% QBOG
Tin_bog = -32 + 273.15; % Inlet BOG temperature [K]
m_bog=610; % Gas mass flow [Kg/h] (Neskaa, 2010)
[ QHA_07, n_bog] = BOGC( Tin_bog, m_bog );

%% HA-07 Simplified
% Global energy balance
% This energy is adjusted to the requirements of the natural gas
% by modifying the refrigerant flow. These requirements are in BOG.m.
z(14,:)=z(11,:);
n(14)=n(11);
H(14) = H(11) - QHA_07;
h(14) = H(14)/n(14);

% Expected value
vf(14) = 0.5;

[x(14,:), y(14,:), vf(14), T(14), Zl(14), Zv(14),...
 sl(14), sv(14), beta(14), h(14)] =...
 flashCalEO(z(14,:), P(14), x(12,1:4), y(12,1:4), vf(14), T(13) , Zl(13), Zv(13), ...
            sl(13), sv(13), beta(13) , rho, 'PH', h(14));

%% Mixer
n(15)=n(7)+n(13)+n(14);
z(15,:)= (n(7).*z(7,:)+n(13).*z(13,:)+n(14).*z(14,:))/n(15);
H(15)=H(7)+H(13)+H(14);
h(15)=H(15)/n(15);

% Expected value
vf(15) = 1;

[x(15,:), y(15,:), vf(15), T(15), Zl(15), Zv(15), ...
 sl(15), sv(15), beta(15), h(15)] = ...
 flashCalEO(z(15,:), P(15), x(3,1:4), y(3,1:4), vf(15), T(6) , Zl(3), Zv(3),...
            sl(3), sv(3), beta(3) , rho, 'PH', h(15));

%% QHA_02 second side
z(16,:)=z(15,:);
n(16)=n(15);
H(16)=H(15)+QHA_02;
h(16)=H(16)/n(16);
% Expected value
vf(16) = 1;

[x(16,:), y(16,:), vf(16), T(16), Zl(16), Zv(16),...
 sl(16), sv(16), beta(16), h(16)] = ...
 flashCalEO(z(16,:), P(16), x(15,1:4), y(15,1:4), vf(16), T(17), Zl(15), Zv(15), ...
            sl(15), sv(15), beta(15), rho, 'PH', h(16));


%% Q_env
% It is calculated by carrying out an energy balance
Q_env=H(17)-H(16); % This corresponds to the heat from the environment
z(18,:)=z(17,:);
n(18)=n(17);
H(18)=Q_env+H(16);
h(18)=H(18)/n(18);
h(18)= h(17);

% Expected value
vf(18)= 1;

[x(18,:), y(18,:), vf(18), T(18), Zl(18), Zv(18),...
 sl(18), sv(18), beta(18), h(18)] = ...
 flashCalEO(z(18,:), P(18), x(16,1:4), y(16,1:4), vf(18), T(17), Zl(16), Zv(16), ...
            sl(16), sv(16), beta(16) , rho, 'PH', h(17));
H(18)=n(18)*h(18);

Q_env= H(18)-H(16);

%% Temperature Crosses
Cross=[];
%BOG Temperatures from BOG.m
TBOGin = -32 + 273.15; %Inlet temperature [K]
TBOGout = -154 +273.15; %Outlet temperature [K]
% % HA - 02
Cross(1) = -T(16) + T(5);
Cross(2) = -T(15) + T(6);
% % HA - 04
Cross(3) = - T(13) + T(6);
Cross(4) =  -T(12) + T(8);

% % HA - 07
Cross(5) = - T(14) + TBOGin;
Cross(6) = - T(11) + TBOGout;

%% Output information
% Result = struct('z',[],'x',[],'y',[],'P',[],'T',[],'VF',[],'h',[],'n',[]);

Result.x=x;
Result.y=y;
Result.VF=vf;
Result.P=P;
Result.T=T
Result.h=h;
Result.z=z;
Result.n=n;
Result.Zl = Zl;
Result.Zv = Zv;
Result.sl = sl;
Result.sv = sv;
Result.beta = beta;
% % % Save information to file
% filename='SSresults.mat';
% save(filename)
