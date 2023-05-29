function [D] = hs_(s,X)
  v = s * X;
  D = tanh(v);
 end