#ifndef _QUATERNIONS_H_
#define _QUATERNIONS_H_

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "libmat.h"
#include "structs.h"

//Declaração de Estrutura de Quaternion
typedef struct Quaternion
{float x,y,z,w;} Quaternion;


Quaternion *NewQuarternion(float x, float y, float z, float w);
void quaternion_mul(Quaternion *ql, Quaternion *qr, Quaternion *Qout);
void quaternion_inv(Quaternion *ql, Quaternion *Qout);
void quat2dcm(Quaternion *Q, struct matrix* dcm);
float quatnormalize(Quaternion *Q);
void quatDisplay(Quaternion *Q);







#endif /*_QUATERNIONS_H_*/