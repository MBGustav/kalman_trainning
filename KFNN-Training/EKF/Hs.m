function [out] = Hs(s,X)
%Sem bias!
sum = 1 - tanh(s(1)*X(1) + s(2)*X(2)+s(3)*X(3))^2;

out(1,1) = X(1)*sum;
out(1,2) = X(2)*sum;
out(1,3) = X(3)*sum;
end