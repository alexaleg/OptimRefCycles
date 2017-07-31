function [c ceq] = f(w, Plow, flow)
  %% HA-04 test
  % TODO check model, maybe the heat exchanger is underspecified.
  % Parameters
  NC = 5;
  z = [0.1906    0.5145    0.2700    0.0242    0.0007];
  Phigh = 20e5;

  % Stream 6
  T6 = 198.15;
  h6 = -8.1426e3;
  H6 = h6*flow;

  % Extract data
  stream = zeros(3,3*NC);
  for i = 0:2
    stream((i+1),:) = w(i*3*NC+1:(i+1)*3*NC);
  end
  VF8 = 0;
  % Sf = w(end);
  Q = w(end)
  % Sf = 0.65;
  Sf = w(end-1);

  % Stream 8
  w8 = stream(1,:);
  h8 = h6 - Q/flow;
  [c1, ceq1, hf8, hg8]=PVCal(w8, z, Phigh, VF8);
  h8_c = hf8*(1-VF8)+hg8*VF8;
  ceq8 = h8 - h8_c;

  % Stream 12
  VF12 = 0;
  w12 = stream(2,:);
  [c2, ceq2] = PHCal(w12, z, Plow, h8);
  flow12 = flow*Sf;
  H12 = flow12*h8;

  % Stream 13
  w13 = stream(3,:);
  h13 = h8 + Q/flow12;
  [c3, ceq3]=PHCal(w13, z, Plow, h13);

  % Add a condition to optimize the split flow
  cSf = w13(9)-1;

  % Constraints
  c = [c1; c2; c3];
  ceq = [ceq1; ceq8; ceq2; ceq3; cSf];
