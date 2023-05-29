#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>
#include <assert.h>
#include <stdarg.h>
#include "errors.h"
#include "util.h"
#include "structs.h"
#include "linear_algebra.h"
#include "linear_algebra_object.h"

struct vector *vector_new(int length)
{
    assert(length >= 0);

    struct vector *new_vector = malloc(sizeof(struct vector));
    check_memory((void *)new_vector);

    DATA(new_vector) = malloc((sizeof(double)) * length);
    check_memory((void *)DATA(new_vector));

    new_vector->length = length;
    OWNS_MEMORY(new_vector) = true;
    MEMORY_OWNER(new_vector) = NULL;
    REF_COUNT(new_vector) = 0;

    return new_vector;
}

// Copy all the data in a given vector into a new vector
void vector_copy_into(struct vector *reciever, struct vector *v)
{
    assert(v->length == reciever->length);
    for (int i = 0; i < v->length; i++)
    {
        VECTOR_IDX_INTO(reciever, i) = VECTOR_IDX_INTO(v, i);
    }
}

struct matrix *matrix_new(int n_row, int n_col)
{
    assert(n_row >= 1 && n_col >= 1);
    struct matrix *new_matrix = malloc(sizeof(struct matrix));
    check_memory((void *)new_matrix);

    DATA(new_matrix) = malloc((sizeof(double)) * n_row * n_col);
    check_memory((void *)DATA(new_matrix));

    new_matrix->n_row = n_row;
    new_matrix->n_col = n_col;
    OWNS_MEMORY(new_matrix) = true;
    // MEMORY_OWNER(new_matrix) = NULL;
    REF_COUNT(new_matrix) = 0;

    return new_matrix;
}

struct vector *vector_new_view(struct linalg_obj *parent, double *view, int length)
{
    assert(length >= 0);
    // TODO: Make this view check work.
    /* Check that pointers to the beginning and end of view vector live
       within the data segment of the parent object.

       This doesn't work because matricies have no length.  This could be
       a property of linalg_obj, but then a macro would be needed to
       make the lookup type generic.

    assert(DATA(parent) <= view && view < DATA(parent) + parent->length);
    assert(view + length <= DATA(parent) + parent->length);
    */

    struct vector *new_vector = malloc(sizeof(struct vector));
    check_memory((void *)new_vector);

    DATA(new_vector) = view;
    new_vector->length = length;
    OWNS_MEMORY(new_vector) = false;
    MEMORY_OWNER(new_vector) = parent;
    REF_COUNT(new_vector) = 0;
    REF_COUNT(parent) += 1;

    return new_vector;
}

/* Print a vector to the console like:
    [1, 2, 3, 4]
*/
// TODO: Maybe this should return a string, we are computing the representation and
// displaying it in the same method.
void vector_print(struct vector *v)
{

    if (v->length == 0)
    {
        printf("[]\n");
    }
    else if (v->length == 1)
    {
        printf("[%.2lf]\n", VECTOR_IDX_INTO(v, 0));
    }
    else
    {
        printf("[%.2lf", VECTOR_IDX_INTO(v, 0));
        for (int i = 1; i < v->length - 1; i++)
        {
            printf(", ");
            printf("%.2lf", VECTOR_IDX_INTO(v, i));
        }
        printf(", %.2lf]\n", VECTOR_IDX_INTO(v, v->length - 1));
    }
}

struct vector *matrix_row_view(struct matrix *M, int row)
{
    assert(0 <= row && row <= M->n_row - 1);
    // double *row_p = DATA(M) + (row * M->n_col);
    // struct vector *r = vector_new_view((struct linalg_obj *)M, row_p, M->n_col);
    struct vector *r = vector_new(M->n_col);
    for (int i = 0; i < M->n_col; i++)
    {
        VECTOR_IDX_INTO(r, i) = MATRIX_IDX_INTO(M, row, i);
    }
    return r;
}
struct vector *matrix_col_view(struct matrix *M, int col)
{
    // boundary check
    assert(0 <= col && col < M->n_col);

    struct vector *r = vector_new(M->n_row);

    for (int i = 0; i < M->n_col; i++)
        VECTOR_IDX_INTO(r, i) = MATRIX_IDX_INTO(M, i, col);
    return r;
}
void vector_free(struct vector *v)
{
    struct linalg_obj *mem_owner;
    if (OWNS_MEMORY(v))
    {
        if (REF_COUNT(v) == 0)
        {
            free(DATA(v));
            free(v);
        }
        else
        {
            raise_non_zero_reference_free_error();
        }
    }
    else
    {
        if (REF_COUNT(v) == 0)
        {
            mem_owner = MEMORY_OWNER(v);
            REF_COUNT(mem_owner) -= 1;
            free(v);
        }
        else
        {
            raise_non_zero_reference_free_error();
        }
    }
}

void matrix_free(struct matrix *v)
{
    struct linalg_obj *mem_owner;
    if (OWNS_MEMORY(v))
    {
        if (REF_COUNT(v) == 0)
        {
            free(DATA(v));
            free(v);
        }
        else
        {
            raise_non_zero_reference_free_error();
        }
    }
    else
    {
        if (REF_COUNT(v) == 0)
        {
            mem_owner = MEMORY_OWNER(v);
            REF_COUNT(mem_owner) -= 1;
            free(v);
        }
        else
        {
            raise_non_zero_reference_free_error();
        }
    }
}

/* Print a matrix to the console like:
    [
      [1, 2],
      [3, 4]
    ]
*/

void matrix_print(struct matrix *M)
{
    struct vector *current_row;
    printf("[\n");
    for (int i = 0; i < M->n_row; i++)
    {
        printf("  ");
        current_row = matrix_row_view(M, i);
        vector_print(current_row);
        vector_free(current_row);
    }
    printf("]\n");
}

// Create a new matrix of a given dimension filled with zeros
struct matrix *matrix_zeros(int n_row, int n_col)
{
    assert(n_row >= 1 && n_col >= 1);
    struct matrix *M = matrix_new(n_row, n_col);
    for (int i = 0; i < n_row * n_col; i++)
    {
        DATA(M)
        [i] = 0;
    }
    return M;
}

// Create a new vector of a given dimension filled with zeros
struct vector *vector_zeros(int elem)
{
    assert(elem >= 1);
    struct vector *V = vector_new(elem);

    // memset -> set all positions in 0 (
    memset(DATA(V), 0, elem * sizeof(double));

    return V;
}

// Create a (square) identity matrix of a given size
struct matrix *matrix_identity(int size)
{
    assert(size >= 1);
    struct matrix *M = matrix_new(size, size);
    for (int i = 0; i < size * size; i++)
    {
        if (MATRIX_ROW(M, i) == MATRIX_COL(M, i))
        {
            DATA(M)
            [i] = 1;
        }
        else
        {
            DATA(M)
            [i] = 0;
        }
    }
    return M;
}

/* Transpose a matrix.
   Note that this creates a new matrix, and copies the data into the new matrix.
*/
struct matrix *matrix_transpose(struct matrix *M)
{
    struct matrix *Mt = matrix_new(M->n_col, M->n_row);
    for (int i = 0; i < M->n_row; i++)
    {
        for (int j = 0; j < M->n_col; j++)
        {
            MATRIX_IDX_INTO(Mt, j, i) = MATRIX_IDX_INTO(M, i, j);
        }
    }
    return Mt;
}

// Compute the matrix product of two aligned matricies
struct matrix *matrix_multiply(struct matrix *Mleft, struct matrix *Mright)
{
    assert(Mleft->n_col == Mright->n_row);
    struct matrix *Mprod = matrix_zeros(Mleft->n_row, Mright->n_col);
    for (int i = 0; i < Mprod->n_row; i++)
    {
        for (int k = 0; k < Mleft->n_col; k++)
        {
            for (int j = 0; j < Mprod->n_col; j++)
            {
                MATRIX_IDX_INTO(Mprod, i, j) +=
                    MATRIX_IDX_INTO(Mleft, i, k) * MATRIX_IDX_INTO(Mright, k, j);
            }
        }
    }
    return Mprod;
}

void matrix_multiply_into(struct matrix *reciever, struct matrix *Mleft, struct matrix *Mright)
{
    assert(Mleft->n_col == Mright->n_row);
    // Zero out the reciever matrix.
    for (int i = 0; i < reciever->n_row; i++)
    {
        for (int j = 0; j < reciever->n_col; j++)
        {
            MATRIX_IDX_INTO(reciever, i, j) = 0;
        }
    }
    // Now multiply.
    for (int i = 0; i < Mleft->n_row; i++)
    {
        for (int k = 0; k < Mleft->n_col; k++)
        {
            for (int j = 0; j < Mright->n_col; j++)
            {
                MATRIX_IDX_INTO(reciever, i, j) +=
                    MATRIX_IDX_INTO(Mleft, i, k) * MATRIX_IDX_INTO(Mright, k, j);
            }
        }
    }
}

// Compute the product of an aligned matrix vector pair
struct vector *matrix_vector_multiply(struct matrix *M, struct vector *v)
{
    assert(M->n_col == v->length);
    struct vector *w = vector_new(M->n_row);
    // w->data = (double*) malloc(sizeof(double) * w->length);
    double sum;
    for (int i = 0; i < M->n_row; i++)
    {
        sum = 0.0;
        for (int j = 0; j < M->n_col; j++)
            sum += MATRIX_IDX_INTO(M, i, j) * VECTOR_IDX_INTO(v, j);
        VECTOR_IDX_INTO(w, i) = sum;
    }

    return w;
}

// Create a (square) diagonal matrix of a given size
struct matrix *matrix_diagonal(struct vector *v, int size)
{
    int cont = 0;
    assert(size >= 1);
    struct matrix *M = matrix_new(size, size);
    for (int i = 0; i < size * size; i++)
    {
        if (MATRIX_ROW(M, i) == MATRIX_COL(M, i))
        {
            DATA(M)
            [i] = VECTOR_IDX_INTO(v, cont);
            cont++;
        }
        else
        {
            DATA(M)
            [i] = 0;
        }
    }
    return M;
}

// struct matrix* cholesky(struct matrix *A, int scal){
// /*Matrix Decomposition: Cholesky
//     Input: A(must NxN)
//            scal(value to multiply before inversion)

//     Output: L(Inverted Matrix) & return null case not invertible
// */
//     //so eh possivel inversao de matrizes quadradas
//     assert(A->n_col == A->n_row);
//     int n = A->n_col;

//     struct matrix* L = matrix_new(n, n);
//     if (L == NULL)
//         exit(EXIT_FAILURE);

//     for (int i = 0; i < n; i++)
//         for (int j = 0; j < (i+1); j++) {
//             double s = 0;
//             for (int k = 0; k < j; k++)
//                 // s += L[i * n + k] * L[j * n + k];
//                 s += MATRIX_IDX_INTO(L,i,k) * MATRIX_IDX_INTO(L,j,k);

//             // L[i * n + j] = (i == j) ?
//             //                sqrt(A[i * n + i] - s) :
//             //                (1.0 / L[j * n + j] * (A[i * n + j] - s));
//             if(i == j) MATRIX_IDX_INTO(L,i,j) = sqrt(MATRIX_IDX_INTO(A,i,i) - s);
//             else MATRIX_IDX_INTO(L,i,j) = (1.0 / MATRIX_IDX_INTO(L, j, j) * (MATRIX_IDX_INTO(A, i, j) - s));
//         }

//     return L;
// }

struct matrix *cholesky(struct matrix *A, double scal)
{
    int i, j, k;
    double sum;
    struct matrix *G;
    double var;
    int n = A->n_col;
    /* CREATE THE MATRIX WITH THE SAME SIZE OF A */
    G = matrix_new(n, n);

    for (k = 0; k < n; k++)
        for (i = 0; i <= k; i++)
        {
            sum = 0;
            for (j = 0; j < i; j++)
            // sum += G[i][j] * G[k][j];
            {
                sum += MATRIX_IDX_INTO(G, i, j) * MATRIX_IDX_INTO(G, k, j);
            }
            if (i == j)
            {
                assert(scal * MATRIX_IDX_INTO(A, i, i) - sum > 0);
                // printf("sum: %.2lf, %.2lf, %.2lf", scal, MATRIX_IDX_INTO(A, i, i), sum);
                var = sqrt((scal * MATRIX_IDX_INTO(A, i, i)) - sum);
                MATRIX_IDX_INTO(G, i, i) = var;
            }
            else
            {
                var = (1.0 / MATRIX_IDX_INTO(G, i, i)) * (scal * MATRIX_IDX_INTO(A, k, i) - sum);
                MATRIX_IDX_INTO(G, k, i) = var;
            }
        }

    // matrix_free(A);
    return G;
}

void matrix_scal_multiply(struct matrix *A_in, struct matrix *B_out, double scal)
{
    assert(A_in != NULL && B_out != NULL);
    // TODO: should i verify width && height?

    for (int i = 0; i < A_in->n_row; i++)
        for (int j = 0; j < A_in->n_col; j++)
            MATRIX_IDX_INTO(B_out, i, j) = scal * MATRIX_IDX_INTO(A_in, i, j);
}

struct matrix *vector_vector_multiply(struct vector *A, struct vector *B)
{

    // check memory
    assert(A != NULL && B != NULL);

    // Math Operation: A[M,1] * B[1,N] = C[M,N]
    struct matrix *out = matrix_zeros(A->length, B->length);

    // In this function, we're supposing B already transposed..
    int m = A->length;
    int n = B->length;

    for (int i = 0; i < m; i++)
        for (int j = 0; j < n; j++)
            MATRIX_IDX_INTO(out, i, j) = VECTOR_IDX_INTO(A, i) * VECTOR_IDX_INTO(B, j);

    return out;
}

void vector_sum(struct vector *v1, struct vector *v2, struct vector *vout)
{

    // check if they're not null...
    assert(v1 != NULL && v2 != NULL && vout != NULL);

    // check if they have same length
    assert(v1->length == v2->length && v2->length == vout->length);

    // time to sum
    for (int i = 0; i < v1->length; i++)
        VECTOR_IDX_INTO(vout, i) = VECTOR_IDX_INTO(v1, i) + VECTOR_IDX_INTO(v2, i);
}

void vector_scalar_multiply(struct vector *vin, struct vector *vout, double alpha)
{
    assert(vin->length == vout->length);

    for (int i = 0; i < vin->length; i++)
        VECTOR_IDX_INTO(vout, i) = alpha * VECTOR_IDX_INTO(vin, i);
}

void matrix_sum(struct matrix *Ml, struct matrix *Mr, struct matrix *Mout)
{

    // check memory
    assert(Ml != NULL && Mr != NULL && Mout != NULL);

    // check size
    assert(MATRIX_SIZE(Ml) == MATRIX_SIZE(Mr) && MATRIX_SIZE(Mout));

    for (int i = 0; i < Ml->n_col; i++)
        for (int j = 0; j < Ml->n_row; j++)
            MATRIX_IDX_INTO(Mout, i, j) = MATRIX_IDX_INTO(Ml, i, j) + MATRIX_IDX_INTO(Mr, i, j);
}

void matrix_vector_ColSum(struct matrix *M, struct vector *V, int col)
{

    // Boundary Checking
    assert(0 <= col && col < M->n_col && V->length == col);
    // Memory Check
    check_memory(DATA(M));
    check_memory(DATA(V));

    for (int i = 0; i < M->n_row; i++)
        MATRIX_IDX_INTO(M, i, col) = VECTOR_IDX_INTO(V, i);
}
void matrix_vector_RowSum(struct matrix *M, struct vector *V, int row)
{

    // mem. check
    assert(M != NULL && V != NULL);
    // Boundary Checking
    assert(0 <= row && row < M->n_row);
    for (int i = 0; i < M->n_col; i++)
        MATRIX_IDX_INTO(M, row, i) = VECTOR_IDX_INTO(V, i);
}

struct matrix* matrix_inverse(struct matrix *A_in)
{

    assert(A_in->n_col == A_in->n_row); // equal sizes --> A
    // assert(B->n_col == B->n_row); // equal sizes --> B
    // assert(B->n_col == A->n_row); // equal sizes --> B && A


    int i, j, k, size = A_in->n_col;
    double temp = 0;

    struct matrix *A = matrix_new(A_in->n_row, A_in->n_col);
    matrix_copy(A_in, A);
    struct matrix *I = matrix_identity(size);
    for (k = 0; k < size; k++){
        temp = MATRIX_IDX_INTO(A,k, k);
        for (j = 0; j < size; j++){
            MATRIX_IDX_INTO(A,k, j) /= temp;
            MATRIX_IDX_INTO(I, k, j)/= temp;
        }
        for (i = 0; i < size; i++){
            temp = MATRIX_IDX_INTO(A,i, k);
            for (j = 0; j < size; j++){
                if (i == k) break;
                MATRIX_IDX_INTO(A, i, j) -= MATRIX_IDX_INTO(A, k, j) * temp;
                MATRIX_IDX_INTO(I, i, j) -= MATRIX_IDX_INTO(I, k, j)* temp;
            }
        }
    }
    matrix_free(A);
    return I;
    
}

void matrix_copy(struct matrix *in, struct matrix * out){
    
    assert(in->n_col == out->n_col);
    assert(in->n_row == out->n_row);

    int col = in->n_col;
    int row = in->n_row;

    for(int i = 0; i < col; i++)
        for(int j = 0; j < row; j++)
            MATRIX_IDX_INTO(out, i,j) = MATRIX_IDX_INTO(in, i,j);
}