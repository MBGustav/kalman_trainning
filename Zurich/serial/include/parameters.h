#ifndef _PARAMETER_H_
#define _PARAMETER_H_

#define d_type double
#define PI 3.14159265358979


#define angle2rad PI/180
#define rad2angle 180/PI
#define tau  100  // 
#define g -9.81


//Definition gyro
#define  gyro_x 1e-7
#define  gyro_y 1e-7
#define  gyro_z 1e-7
#define  gyro_bias_x 0.0000000001
#define  gyro_bias_y 0.0000000001
#define  gyro_bias_z 0.0000000001

//Definition Accelerometer 
#define  acce_x 0.01
#define  acce_y 0.01
#define  acce_z 0.01
#define  acce_bias_x 0.000000001
#define  acce_bias_y 0.000000001
#define  acce_bias_z 0.000000001

//Define Magnetometer
#define  magn_x 0.1
#define  magn_y 0.1
#define  magn_z 0.1

#define  pos_gps_x 0.0001
#define  pos_gps_y 0.0001
#define  pos_gps_z 0.0001
#define  vel_gps_x 1
#define  vel_gps_y 1
#define  vel_gps_z 1 



#endif /*_PARAMETER_H_*/