function [landmarks] =  read_world(world)
    % L� a defini��o do mundo e retorna uma lista de pontos de refer�ncia, 
    % ou seja o 'mapa'
    %
    % O dicion�rio retornado cont�m uma lista de pontos de refer�ncia, cada
    % um com a seguinte informa��o: [id, x, y]

    landmarks(:,1) = world.VarName1(1:end);
    landmarks(:,2) = world.VarName2(1:end);
    landmarks(:,3) = world.VarName3(1:end);
end