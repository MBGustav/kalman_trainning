#ifndef _STRUCTS_H__
#define _STRUCTS_H__

#include <stdbool.h>
#include "linear_algebra_object.h"

// ################################################################
// ########################## Structures ##########################
// ################################################################


// #################### General Purpose Arrays ####################

struct matrix
{
    struct linalg_obj la_obj;
    int n_row;
    int n_col;
};

struct vector {
    struct linalg_obj la_obj;
    int length;
};

/*

// ######################## IMU's structs #########################
struct axes
{
    double x;
    double y;
    double z;
};

struct IMU
{
    struct axes acel;
    struct axes rate;
    struct axes mag;
};

struct data_IMU
{
    struct IMU bin;
    struct IMU scaled;
    struct IMU adj;
};

struct euler
{
    double yaw;
    double pitch;
    double roll;
};

// ##################### State Space's structs ####################
*/

struct wm
{
    struct matrix *Q;
    struct matrix *Qeta;
    struct matrix *Qeta_sr;
    struct matrix *R;
    struct matrix *R_sr;
    struct matrix *P;
    struct matrix *P_sr;
    struct matrix *Pp;
    
};

struct ss_d
{

    struct matrix *Phi;
    struct matrix *G;
    struct matrix *H;
    struct wm wm;

    //sempre vector ou devemos usar como matrix? sempre vector !
    struct vector *xm;
    struct vector *xp;
    // struct vector *z;
    struct vector *x;
};

#endif // _STRUCTS_H__