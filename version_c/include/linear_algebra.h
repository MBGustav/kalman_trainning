#ifndef _LINEAR_ALGEBRA_H__
#define _LINEAR_ALGEBRA_H__

#include <stdarg.h>
#include "linear_algebra_object.h"

#ifndef _MATRIX_MACROS
#define _MATRIX_MACROS
#define MATRIX_ROW(M, i) ((i) / (M->n_row))
#define MATRIX_COL(M, i) ((i) % (M->n_row))
#define MATRIX_SIZE(M)   ((M->n_row * M->n_col))
#define MATRIX_IDX(M, r, c) ((M->n_col * r + c))

#define MATRIX_IDX_INTO(M, r, c) (DATA(M)[MATRIX_IDX(M, r, c)])
#endif

#ifndef _VECTOR_MACROS
#define _VECTOR_MACROS
#define VECTOR_IDX_INTO(v, i) (DATA(v)[i])
#endif

struct vector* vector_new(int length);
struct matrix* matrix_new(int n_row, int n_col);
struct matrix* matrix_zeros(int n_row, int n_col);
struct vector* vector_zeros(int elem);
struct matrix* matrix_identity(int size);
struct matrix* matrix_diagonal(struct vector *v, int size);
struct matrix* matrix_transpose(struct matrix *M);
struct matrix* matrix_multiply(struct matrix* Mleft, struct matrix* Mright);
struct vector* matrix_vector_multiply(struct matrix *M, struct vector *v);
struct vector* vector_new_view(struct linalg_obj *parent, double *view, int length);
void vector_print(struct vector *v);
struct vector *matrix_row_view(struct matrix *M, int row);
void vector_copy_into(struct vector *reciever, struct vector *v);
void matrix_print(struct matrix *);
void matrix_multiply_into(struct matrix *reciever, struct matrix *Mleft, struct matrix *Mright);

void vector_free(struct vector *v);

//TODO:to implement.. OK
void matrix_free(struct matrix *m);

//TODO:to implement.. OK
struct matrix* cholesky(struct matrix* A, double scal);

void matrix_scal_multiply(struct matrix *A_in, struct matrix* B_out, double scal);

//TODO: to implement.. OK 
struct matrix* vector_vector_multiply(struct vector* A, struct vector* B);

//TODO: scalar sum
void vector_sum(struct vector* v1, struct vector* v2, struct vector* vout);

//TODO: collumn view.. OK
struct vector *matrix_col_view(struct matrix *M, int row);


//TODO: scalar vector multiply
void vector_scalar_multiply(struct vector* vin, struct vector* vout, double alpha);


//TODO: matrix sum
void matrix_sum(struct matrix* Ml, struct matrix* Mr, struct matrix *Mout);


void matrix_vector_ColSum(struct matrix *M, struct vector *V, int col);
void matrix_vector_RowSum(struct matrix *M, struct vector *V, int row);

struct matrix * matrix_inverse(struct matrix *A);

void matrix_copy(struct matrix *in, struct matrix * out);


#endif // __LINEAR_ALGEBRA_H__