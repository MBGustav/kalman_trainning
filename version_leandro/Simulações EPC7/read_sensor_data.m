function [sensor_readings]= read_sensor_data(sensordata)
    
    % L� a odometria e as leituras do sensor de um arquivo .mat
    %
    % Os dados s�o retornados em um dicion�rio, onde u_t e z_t s�o 
    % armazenados juntos da seguinte forma:
    % {odometria, sensor}
    %
    % onde "odometria" tem os campos r1, r2, t que cont�m os valores de
    % as vari�veis do modelo de movimento com nome id�ntico, e sensor � uma
    % lista de leituras do sensor com id, faixa e rumo como valores
    %
    % Os valores da odometria e do sensor s�o acessados da seguinte forma:
    % odometry_data = sensor_reading(timestep).odometria
    % sensor_data = sensor_readings(timestep).sensor
    
    % Tamanho do vetor de dados
    N = height(sensordata);
    
    timestep = 1;
    for k = 1:N
        
        if strcmp(string(sensordata.ODOMETRY(k)),'ODOMETRY')
            sensor_readings(timestep).odometry.r1 = sensordata.VarName2(k);
            sensor_readings(timestep).odometry.t = sensordata.VarName3(k);
            sensor_readings(timestep).odometry.r2 = sensordata.VarName4(k);
            timestep = timestep + 1;
            sensorstep = 1;
            
        end
        if strcmp(string(sensordata.ODOMETRY(k)),'SENSOR')
            sensor_readings(timestep-1).sensor.id(sensorstep) = sensordata.VarName2(k);
            sensor_readings(timestep-1).sensor.range(sensorstep) = sensordata.VarName3(k);
            sensor_readings(timestep-1).sensor.bearing(sensorstep) = sensordata.VarName4(k);
            sensorstep = sensorstep+1;       
        end
    end