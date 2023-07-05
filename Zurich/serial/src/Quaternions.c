
#include "Quaternions.h"


Quaternion *NewQuarternion(float x, float y, float z, float w){
    Quaternion *Qout = (Quaternion *)malloc(sizeof(Quaternion));

    if(Qout == NULL)
        return NULL;


    Qout->w = w;
    Qout->x = x;
    Qout->y = y;
    Qout->z = z;

    return Qout;
}

void quaternion_mul(Quaternion *ql, Quaternion *qr, Quaternion *Qout){

    if(ql == NULL || qr == NULL || Qout == NULL)
        fprintf(stderr, "Could not multiply!\n");
    
    Qout->w = ql->w * qr->w - ql->x * qr->x - ql->y * qr->y - ql->z * qr->z;
    Qout->x = ql->w * qr->x + ql->x * qr->w + ql->y * qr->z - ql->z * qr->y;
    Qout->y = ql->w * qr->y + ql->y * qr->w + ql->z * qr->x - ql->x * qr->z;
    Qout->z = ql->w * qr->z + ql->z * qr->w + ql->x * qr->y - ql->y * qr->x;   
}

void quat2dcm(Quaternion *Q, struct matrix* dcm){

    float norm = quatnormalize(Q); 
    float x,y,z,w;

    x = Q->x;/*qn(:,1)*/
    y = Q->y;/*qn(:,2)*/
    z = Q->z;/*qn(:,3)*/
    w = Q->w;/*qn(:,4)*/

    dcm->elements[3*0 + 0] = x*x + y*y - z*z - w*w;
    dcm->elements[3*0 + 1] = 2*(y*z + x*w);
    dcm->elements[3*0 + 2] = 2*(y*w - x*z);
    dcm->elements[3*1 + 0] = 2*(y*z - x*w);
    dcm->elements[3*1 + 1] = x*x - y*y + z*z - w*w;
    dcm->elements[3*1 + 2] = 2*(z*w + x*y);
    dcm->elements[3*2 + 0] = 2*(y*w + x*z);
    dcm->elements[3*2 + 1] = 2*(z*w - x*y);
    dcm->elements[3*2 + 2] = x*x - y*y - z*z + w*w;
}



void quaternion_inv(Quaternion *ql, Quaternion *Qout){

    if(ql == NULL || Qout == NULL){
        fprintf(stderr, "Could not multiply!\n");
        exit(1);
    }
    float q0 = ql->w;
    float q1 = ql->x;    
    float q2 = ql->y;    
    float q3 = ql->z;    

    float norm_q = sqrt(q0 * q0 + q1 * q1 + q2 * q2 + q3 * q3);
    Qout->w =  q0 / norm_q;
    Qout->x = -q1 / norm_q;
    Qout->y = -q2 / norm_q;
    Qout->z = -q3 / norm_q;

}


float quatnormalize(Quaternion *Q){
    if(Q == NULL){
        fprintf(stderr, "Could not QuatNormalize!\n");
        exit(1);
    }
    float norm = Q->x*Q->x + Q->y*Q->y + Q->z*Q->z + Q->w*Q->w;
    return norm;
}

void quatDisplay(Quaternion *Q){
    if(Q == NULL ){
        fprintf(stderr, "Could not Display!\n");
        exit(1);
    }
    printf("\nQ = %.3f, %.3f, %.3f, %.3f \n", Q->x, Q->y, Q->z, Q->w);
}
