function [J]=f(w, rho)
  VF = w(9);
  sl = w(12);
  sv = w(13);
  J = penalty(VF, sl, sl, rho);
