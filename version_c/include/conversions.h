
#ifndef __CONVERSIONS_H__
#define __CONVERSIONS_H__


void angle2dcm(struct matrix *ptr_dcm, struct vector *ptr_angles);
void enu2ned(struct vector *ptr_axis);
// void geodetic2enu(struct vector *ptr_ENU, struct vector *ptr_curr_geodetic, struct vector *ptr_init_geodetic);
void geodetic2ned(struct vector *ptr_ENU, struct vector *ptr_curr_geodetic, struct vector *ptr_init_geodetic);


#endif //__CONVERSIONS_H__