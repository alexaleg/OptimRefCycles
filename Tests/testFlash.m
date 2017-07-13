clear all
z = [0.025 0.65 0.15 0.15 0.025];
z = [0.1074 0.3292 0.4096 0.1345 0.0193];
P = 55e5;
P = 3e5
sl = 0;
sv = 0;
beta = 1;
Zl = 0.25;
Zv = 0.65;
T = 200;
x0 = z(1:end-1); y0=x0; VF0 = 0.74;
rho = 5e4;
listT = [ 300 296 292.882 270 245 212.8 209.5 206.5];
listVF = zeros(size(listT));
listH = listVF;

% Feed initialization
sprintf('PT')
[x0, y0, VFt, Tt, Zl0, Zv0, sl0, sv0, beta0, ht] = flashCalEO(z, P, x0, y0, VF0, T, Zl, Zv, sl, sv, beta, rho, 'PT');
sprintf('PH')
[x0, y0, VFh, Th, Zl0, Zv0, sl0, sv0, beta0, hh] = flashCalEO(z, P, z(1:4), z(1:4), VF0, T, Zl, Zv, sl, sv, beta, rho, 'PH', ht);
VFt
VFh
Tt
Th
% listVF(1)=VF;
% listH(1)=h;
% for i=2:length(listT)
%   T = listT(i )
%   [x, y, VF, T, Zl, Zv, sl, sv, beta,h] = flashCalEO(z, P, x0(1:4), y0(1:4), VF0, T, Zl0, Zv0, sl0, sv0, beta0, rho, 'PT');
%   VF
%   listVF(i)=VF;
%   listH(i)=h;
% end
% listVF
% listTH = zeros(size(listT));
% err = listTH;
% listVFH = listTH;
% for i=1:length(listT)
%   [x0, y0, VF, T, Zl0, Zv0, sl0, sv0, beta0, h] = flashCalEO(z, P, x0(1:4), y0(1:4), VF0, T, Zl, Zv, sl, sv, beta, rho, 'PH', listH(i));
%   listTH(i)=T;
%   err(i)= listH(i)-h;
%   listVFH(i)=VF;
% end
% listT
% listTH
