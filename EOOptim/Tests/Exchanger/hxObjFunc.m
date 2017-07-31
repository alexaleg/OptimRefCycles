function [J] = f(w,rho)
  % Extract data
  NC = 5;
  VF = zeros(1,3);
  sl = zeros(1,3);
  sv = zeros(1,3);

  for i=1:3
    VF(i) = w((2*NC-1)+(i-1)*3*NC);
    sl(i) = w((2*NC+2)+(i-1)*3*NC);
    sv(i) = w((2*NC+3)+(i-1)*3*NC);
  end
  VF(1);
  % J8  = (VF(1)+sl(1))*rho;
  J8 = penalty(VF(1),sl(1),sv(1),rho);
  J12 = penalty(VF(2),sl(2),sv(2),rho);
  J13 = penalty(VF(3),sl(3),sv(3),rho);

  J = J8 + J12 + J13;
  % J = J8 + J12;
