
x = [0.1074 0.3292 0.4096 0.1345 0.0193];
T1 = -42 + 273.15;
P = 3e5;
phase = 2;
[Z,phi,H,V,S]= srk(x,T1,P,phase);

T2 = 67.25 + 273.15;
P2 = 20e5;

T0 = 390;
Tout = fzero(@(T)Sentropy( S, x ,T, P2 ),T0)


[Z,phig,H2,V,S2]= srk(x,Tout,P2,phase);

[Z,phil,H2l,V,S2l]= srk(x,Tout,P2,1);
eta = 0.74;
W = (H2 - H)/eta

H2_real = H + W;

T3 = 90.3 + 273.15;
P3 = 20e5;

Tout1 = fzero(@(T)Centhalpy( H2_real, x ,T, P2 ),T0)
[Z,phi,H3,V,S3]= srk(x,T3,P3,phase);
