function [c,ceq]=f(w, Pin, Pout)
    %% Extract data
  stream = w(1:15);
  y = [w(5:8) 1-sum(w(5:8))];

  Tout = w(end);
  Tout_real = w(end-1);

  % Iteration on Tout
  z = [0.1074 0.3292 0.4096 0.1345 0.0193];
  Tin = 236.15;

  %% Curves calculation
  p = [0.0015   -0.0314    0.1911    0.4217]; %3rd degree polynomial.
  rel = Pout/Pin; % Values from 2 to 10
  eta = polyval(p,rel);

  [~,~,Hi,~,Si]= srk(z,Tin,Pin,2);                   % Inlet enthalpy and entropy
  [~,~,Ho,~,So]= srk(z,Tout,Pout,2);                  % Output isoentropic enthalpy
  f1 = (Si-So)*1e3;

  W = (Ho-Hi)/eta;    % Work calculation for
  Ho_real = Hi + W;   % Real output enthalpy


  [c,ceq]=PHCal(w, z, Pout, Ho_real);

  ceq = [ceq; f1];
