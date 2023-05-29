function [D] = hs_(s,X)
  W = s; 
  dot = s(1)*X(1) + s(2)*X(2) + s(3)*X(3);
  D = tanh(dot);
  
 end