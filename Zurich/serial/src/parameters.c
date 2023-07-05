#include <stdio.h>
#include <stdlib.h>
#include "parameters.h"
#include "libmat.h"
#include "structs.h"

void *parameters(struct wm *system){


    if(system == NULL){
        printf("Error in parameters");
        exit(1);
    }

    struct matrix *Q = matrix_m(12, 12);
    struct matrix *R = matrix_m(12, 12);

    eye_m(0, Q);
    eye_m(0, R);
    
    //Ruido do Processo
    
    Q->elements[12*0 +0] = gyro_x;
    Q->elements[12*1 +1] = gyro_y;
    Q->elements[12*2 +2] = gyro_z;
    
    Q->elements[12*3 +3] = gyro_bias_x;
    Q->elements[12*4 +4] = gyro_bias_y;
    Q->elements[12*5 +5] = gyro_bias_z;

    Q->elements[12*6 +6] = acce_x;
    Q->elements[12*7 +7] = acce_y;
    Q->elements[12*8 +8] = acce_z;

    Q->elements[12*9 +9 ] = acce_bias_x;
    Q->elements[12*10+10] = acce_bias_y;
    Q->elements[12*11+11] = acce_bias_z;


    //Ruido das medidas
    R->elements[12*0 +0] = magn_x;
    R->elements[12*1 +1] = magn_y;
    R->elements[12*2 +2] = magn_z;

    R->elements[12*3 +3] = acce_x;
    R->elements[12*4 +4] = acce_y;
    R->elements[12*5 +5] = acce_z;

    R->elements[12*6 +6] = vel_gps_x;
    R->elements[12*7 +7] = vel_gps_y;
    R->elements[12*8 +8] = vel_gps_z;

    R->elements[12*9 +9 ] = pos_gps_x;
    R->elements[12*10+10] = pos_gps_y;
    R->elements[12*11+11] = pos_gps_z;



}