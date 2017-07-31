function [c ceq] = f(w, flow)
  %% Stream 8 test

  % Parameters
  NC = 5;
  z = [0.1906    0.5145    0.2700    0.0242    0.0007];
  Phigh = 20e5;

  % Stream 6
  T6 = 198.15;
  h6 = -8.1426e3;
  H6 = h6*flow;

  % Extract data
  stream = zeros(1,3*NC);
  stream(1,:) = w(1:3*NC);

  Q = w(end);
  Sf = 0.65;

  % Stream 8
  H8 = H6 - Q;
  h8 = H8/flow;
  w8 = stream(1,:);
  VF8 = 0;
  [c1, ceq1, hf8, hg8] = PVCal(w8, z, Phigh, VF8);
  h8_c =  hf8*(1-VF8)+hg8*VF8;
  ceq2 = h8_c-h8;
  % Constraints
  c = [c1; ];
  ceq = [ceq1; ceq2];
