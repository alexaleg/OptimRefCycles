clear all
clc
%% Compression Test
% This script tests the reliability of the compressor model and the entropy
% calculations using SRK

% Case 1 - Mini LNG conditions
disp('1')
Pin = 3e5;
Pout = 20e5;

Tin = 231.15;
z1 = [0.1074 0.3292 0.4096 0.1345 0.0193];

N = 1; %Molar flow [kj/kmol]

[ Tout_real ,W,hout,hin] = compIsoS( Pin,Tin,Pout,z1,N);

Tout_hs = 330.9;
W_hs = 4171;
Error(1,1) = 100*abs(Tout_hs-Tout_real)/Tout_hs;
Error(2,1) = 100*abs(W_hs-W)/W_hs;

% Case 2 - Higher pressures
disp('2')
Pin = 20e5;
Pout = 30e5;

Tin = 300;
z1 = [0.1074 0.3292 0.4096 0.1345 0.0193];

N = 1; %Molar flow [kj/kmol]

[ Tout_real ,W,hout,hin] = compIsoS( Pin,Tin,Pout,z1,N);

Tout_hs2 = 324.3;
W_hs2 = 934.7;
Error(1,2) = 100*abs(Tout_hs2-Tout_real)/Tout_hs2;
Error(2,2) = 100*abs(W_hs2-W)/W_hs2;

% Case 3 - Higher pressures + Heavier composition
disp('3')
Pin = 30e5;
Pout = 40e5;

Tin = 400;
z1 = [0.1 0.2 0.1 0.3 0.3];

N = 1; %Molar flow [kj/kmol]

[ Tout_real ,W,hout,hin] = compIsoS( Pin,Tin,Pout,z1,N);

Tout_hs2 = 416.15;
W_hs2 = 904;
Error(1,3) = 100*abs(Tout_hs2-Tout_real)/Tout_hs2;
Error(2,3) = 100*abs(W_hs2-W)/W_hs2;
