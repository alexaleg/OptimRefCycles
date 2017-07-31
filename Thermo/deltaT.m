function [LMTD] = deltaT(Tin_cold, Tout_cold, Tin_hot, Tout_hot)
  %% Funcion to calculate the logarithmic temperature difference
delta1 = Tin_hot - Tout_cold;
delta2 = Tout_hot - Tin_cold;
LMTD = (delta1 - delta2)/log(delta1/delta2);
