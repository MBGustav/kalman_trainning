#ifndef __KFILTER_H__
#define __KFILTER_H__

//Builtin - Libs
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

//Usr - Libraries
#include "linear_algebra.h"
#include "linear_algebra_object.h"
#include "structs.h"
#include "defines.h"
#include "util.h"


/*
    TODO: devo incluir mais uma struct, ou faço diretivas #define para
    uso do filtro EKF <-> UKF ?
    Tomando como base UKF 
*/

typedef struct KFilter
{   
    //Weight Matrices
    struct wm *w;

    // State Space's Data
    struct ss_d *ss;

    //KF parameters
    struct parameters{
        double dT;
        double kpp;

        //TODO: is it better to define this way?
        int nx;
        int nz;
    }p;

    
}KFilter;


/*Initialize matrices:
    1. Initial Syst. State
    2. Initial Stat. Uncert
    3. Set Parameters
*/
KFilter* KF_Init(struct ss_d* pss,struct wm* pw, double dT, double kpp);
// KFilter* KF_Init(double dT, double kpp, int nx, int ny);

void KF_Update(KFilter *KF, struct vector *z);

void KF_SigmaPoints(KFilter *KF, struct matrix* out_Xi, struct vector* W);

void UT(struct matrix *Xi, struct vector *W, struct matrix* noiseCov, struct vector *xm, struct matrix *xcov);

struct vector* fx(struct vector *vin, double dT);

struct vector * hx(const struct vector *in);

#endif //__KFILTER_H__
