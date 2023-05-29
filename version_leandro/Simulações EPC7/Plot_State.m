 function Plot_State(mu, sigma, landmarks, map_limits)
        
        % Posi��es de refer�ncia
        lx = landmarks(:,2);
        ly = landmarks(:,3);
        
        % M�dia da distribui��o normal como estimativa atual
        estimated_pose = mu;

        % C�lculo da elipse de covari�ncia
        covariance = sigma(1:2,1:2);
        [eigenvecs,eigenvals] = eig(covariance);

        % C�lculo do maior autovalor e autovetor
        [max_ind,r] = find(eigenvals == max(max(eigenvals)));
        max_eigvec = eigenvecs(:,max_ind);
        max_eigval = max(max(eigenvals));
        
        min_ind = 1;
        
        if max_ind == 1
            min_ind = 2;
        end
   
        min_eigval = max(eigenvals(:,min_ind));
        min_eigvec= eigenvecs(:,min_ind);
        
        % Valor qui-quadrado para intervalo de confian�a sigma
        chisquare_scale = 1;
        
        % C�lculo da largura e altura da elipse de confian�a
        width = 1 * sqrt(chisquare_scale*max_eigval);
        height = 1 * sqrt(chisquare_scale*min_eigval);
        angle = atan2(max_eigvec(2),max_eigvec(1));
        
        theta_grid = linspace(0,2*pi);
      
        % Elipse em coordenadas x e y
        ellipse_x_r  = width*cos( theta_grid );
        ellipse_y_r  = height*sin( theta_grid );
        
        % Define uma matriz de rota��o
        R = [ cos(angle) sin(angle); -sin(angle) cos(angle) ];

        % Rota��o da elipse em um �ngulo phi
        r_ellipse = [ellipse_x_r;ellipse_y_r]' * R;
        
        [x,y] = pol2cart(estimated_pose(3),1);
        
        figure(1)
        % https://www.visiondummy.com/2014/04/draw-error-ellipse-representing-covariance-matrix/
        % Desenha a elipse do erro
        plot(r_ellipse(:,1) + estimated_pose(1),r_ellipse(:,2) +  estimated_pose(2),'--', 'LineWidth',1.1)
        hold on;
        title('\it Trajet�ria e Orienta��o do Rob�','FontSize',14)
        % % Tra�a os autovetores
        % quiver(estimated_pose(1), estimated_pose(2), max_eigvec(1)*2*sqrt(chisquare_scale*max_eigval), max_eigvec(2)*2*sqrt(chisquare_scale*max_eigval), '-m', 'LineWidth',2);
        % quiver(estimated_pose(1), estimated_pose(2), min_eigvec(1)*2*sqrt(chisquare_scale*min_eigval), min_eigvec(2)*2*sqrt(chisquare_scale*min_eigval), '-g', 'LineWidth',2);
        plot(estimated_pose(1),estimated_pose(2),'o','markersize',10,'MarkerFaceColor','r')
        quiver(estimated_pose(1), estimated_pose(2), 3*x, 3*y,'-r', 'LineWidth',3)
        hold on;
        plot(lx, ly, 'bo','markersize',10,'MarkerFaceColor','b')
        legend({' Erro',' Posi��o',' Orienta��o', '\it Landmarks'},'FontSize',10,'FontName',' Times New Roman','Location','northwest')
        axis(map_limits)
        hold off
        grid on
        xlabel('\it {Localiza��o x}','FontName','Times New Roman','FontSize',14)
        ylabel('\it {Localiza��o y}','FontName','Times New Roman','FontSize',15,'Rotation',90);
        set(gca,'XTick',[-3:13],'YTick',[-2:12],'Box','on');
        set(gcf,'Position',[300 100 600 500],'PaperPositionMode','auto');
        drawnow;
                
    end