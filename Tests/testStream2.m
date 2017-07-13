clear all
z = [0.1074 0.3292 0.4096 0.1345 0.0193];
P = 20e5;
x0 = z(1:end-1);
y0 = x0;
T = 35+273.15;
Zl = 0.2;
Zv = 0.8;
sl0 = 0;
sv0 = 0;
beta0 = 1;
rho = 5e4;
VF0 = 1;
[xcal, ycal, VF, Tcal, Zl, Zv, ~, ~ , ~, h] = flashCalEO(z, P, x0, y0, VF0, T, Zl, Zv , sl0, sv0, beta0, rho, 'PT');
