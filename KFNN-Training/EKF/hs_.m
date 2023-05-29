function [D] = hs_(s,X)
    W = s;
    D = tanh(X*W);
end