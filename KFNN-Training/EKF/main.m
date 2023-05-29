clear all
close all
clc

X = [0 0 1;
     0 1 1;
     1 0 1;
     1 1 1];
 
D = [0;0;1;1];

%%W =  randn(3,1);
W=[-0.47399; 0.73467; 1.5124]

N =  1000;
%% Initial Conditions
ns = size(W,1);
nd = size(D,1);
s = W;
P = eye(ns);
R = 10*eye(nd);
Q = 0.01*eye(ns);

for k = 1:N
 
    %%Predict
    sp = s;
    Pp = P + Q;
    
    %%Update
    z = D;
    for i =1:nd
        H(i,:) = Hs(s,X(i,:));
        h(i,:) = hs_(s,X(i,:));
    end
      
    error = z - h;
    
    S = H*P*H' + R;
    K  = Pp*H'*S^(-1);
    s = sp+K*error;
    P = (eye(ns) -K*H)*Pp;
    cost(k) = sum(abs(error));
end

plot(cost)

 
