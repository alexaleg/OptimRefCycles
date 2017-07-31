function [J] = f(w,rho)
  % Extract data
  NC = 5;
  VF = zeros(1,1);
  sl = zeros(1,1);
  sv = zeros(1,1);

  for i=1
    VF(i) = w((2*NC-1)+(i-1)*3*NC);
    sl(i) = w((2*NC+2)+(i-1)*3*NC);
    sv(i) = w((2*NC+3)+(i-1)*3*NC);
  end
  VF(1);
  J8  = (VF(1)+sl(1))*rho;
  J = J8;
