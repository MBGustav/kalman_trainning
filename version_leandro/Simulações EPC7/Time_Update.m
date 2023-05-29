function [x_pre, y_pre, fXi, hYi, Sigma_pre] =  Time_Update(sensor, odometry, landmarks, xm, Xi, Wm, Wc, Q)
    
    % Atualiza o valor esperado, ou seja, mu e sigma, de acordo com o 
    % modelo de movimento
    %
    % mu: vetor 3x1 que representa a m�dia (x,y,theta) da distribui��o normal 
    % sigma: matriz de covari�ncia 3x3 da distribui��o normal
    
    % Leitura das informa��es de odometria
    delta_rot1 = odometry.r1;
    delta_trans = odometry.t;
    delta_rot2 = odometry.r2;
    delta = [delta_rot1; delta_trans; delta_rot2];
    
    % Propaga��o dos Sigma Points em f[k]
    L = numel(xm);
    fXi = zeros(L,L+1);
    for j=1:2*L+1
        fXi(:,j) = fX(Xi(:,j),delta);
    end
    
    % Pondera��o de Xi
    [x_pre, Sigma_pre] = UTX(fXi, Wm, Wc, Q);
    
    % Propaga��o dos Sigma Points em h[k]
    ids = sensor.id;
    hYi = zeros(numel(ids),L+1);
    for j=1:2*L+1
        hYi(:,j) = hX(Xi(:,j),landmarks,ids);
    end
    
    % Pondera��o de Yi
    [y_pre] = UTY(hYi, Wm);
    
end


%% Fun��es Auxiliares 

% Fun��o de propaga��o dos Sigma Points por f[k]
function Xp = fX(Xi,delta)
    xp(1) = Xi(1) + delta(2)*cos(Xi(3) + delta(1));
    xp(2) = Xi(2) + delta(2)*sin(Xi(3) + delta(1));
    xp(3) = Xi(3) + delta(1) + delta(3);
    
    Xp = [xp(1); xp(2); xp(3)];
end


% Fun��o de propaga��o dos Sigma Points por h[k]
function Yp = hX(Xi,landmarks,ids)

    % Identifica��o do sensor
    tam = numel(ids);
    Yp = zeros(tam,1);
    
    % C�lculo da medi��o de alcance esperada
    for i=1:tam
        for j=1:9
            if ids(i)==j
                lx = landmarks(j,2);
                ly = landmarks(j,3);
                Yp(i,1) = sqrt((lx - Xi(1))^2 + (ly - Xi(2))^2);
            end
        end
    end
    
end


% Fun��o de pondera��o de Xi
function [xp,Pp] = UTX(fXi, Wm, Wc, Q)  

    [n, kmax] = size(fXi);

    xp = zeros(n,1);
    for k=1:kmax
        xp = xp + Wm(k)*fXi(:, k);
    end

    Pp = zeros(n,n);
    for k=1:kmax
        Pp = Pp + Wc(k)*(fXi(:, k) - xp)*(fXi(:, k) - xp)';
    end
    
    Pp = Pp + Q;

end


% Fun��o de pondera��o de Yi
function [yp] = UTY(hYi, Wm)  

    [n, kmax] = size(hYi);

    yp = zeros(n,1);
    for k=1:kmax
        yp = yp + Wm(k)*hYi(:, k);
    end

end