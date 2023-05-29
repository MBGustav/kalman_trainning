clear all 
close all 
clc;


X = [0 0 1; 
     0 1 1;
     1 0 1;
     1 1 1];

D = [0;0;1;1];

%%W =  randn(3,1);
W = [-0.47399; 0.73467; 1.5124];

%%Condicoes iniciais
N = 1;
ns = size(W,1);
nd = size(D,1);
s = W;
kappa = 0;
kmax = 2*ns+1;

%% Modelagem do Sistema
P = eye(ns);
R = eye(nd);
Q = eye(ns);
z = D;

for j=1:N
    %% Measure SigmaPoints
    [Si Ws] = SigmaPoints(s, P, kappa, kmax);
    %%  Estado de Update
     % Propagacao em fs(): Descrição do Sistema
    fSi = zeros(ns, kmax);
    for k=1:kmax
      fSi(:,k) = fs_(Si(:,k));
    end
    
    % Propagacao em hs(): Predição de Saída
    hSi = zeros(nd, kmax);
    for k=1:kmax
      for i =1:nd
        h(i,:) = hs_(Si(:,k), X(i,:));
      end
      hSi(:,k) = h;
    end
    
      % Transformação Uncented -- UT
    [sp Pp] = UT(fSi, Ws, Q);
    [zp Pz] = UT(hSi, Ws, R);
    
    Psz = zeros(ns, nd);
    for k=1:kmax
      % Matriz de Covariancia => Psz
       Psz = Psz + Ws(k) * (fSi(:,k) - sp) * (hSi(:,k) - zp)';
    end
    
    error = z - zp;
    % Estado de Correcao
    K = Psz* Pz^(-1);
    s = sp + K*error;
    P = Pp - K*Pz*K';
    
    cost(j) = sum(abs(error))
end

plot(cost)