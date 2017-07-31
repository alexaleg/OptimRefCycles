function [UA] = f(Q,Tin_cold, Tout_cold, Tin_hot, Tout_hot)
  % Calculate UA
  LMTD = deltaT(Tin_cold, Tout_cold, Tin_hot, Tout_hot);
  UA = Q/LMTD;
