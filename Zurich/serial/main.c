#include <stdio.h>
#include <stdlib.h> 
#include <math.h> 

#include "structs.h"
#include "Quaternions.h"
#include "parameters.h"
#include "libmat.h"


struct matrix * diag(struct matrix **list_M, int total_m){
    
    //acc row (considering square matrix)
    int g_sz = 0;
    for(int i = 0; i < total_m; i++){
        g_sz += list_M[i]->n_col;
    }

    struct matrix *output = matrix_m(g_sz, g_sz);
    eye_m(0.0f, output);
    
    struct matrix *M; 
    int sz;
    int pos = 0;
    for(int idx = 0; idx < total_m; idx++){//for each matrix
        M = list_M[idx];
        sz = M->n_row;
        for(int p_x = 0; p_x<sz ; p_x++){//copy all matrix 
            for(int p_y=0; p_y<sz ; p_y++){
                output->elements[(p_y+pos)*g_sz+(pos+p_x)] =
                M->elements[p_y*sz + p_x];
            }
        }
        pos = sz;
    }
    
    return output;
}

int main(){

    float a = sqrt(4); 

    // Quaternion *Q = NewQuarternion(0,0,0,0);
    // struct matrix *A  = matrix_m(2,2);
    // struct matrix *B  = matrix_m(2,2);
    // A->elements[0] = 1;A->elements[1] = 1;
    // A->elements[2] = 1;A->elements[3] = 1;

    // B->elements[0] = 4;B->elements[1] = 4;
    // B->elements[2] = 4;B->elements[3] = 4;


    // struct matrix *list[2];
    // list[0] = A;
    // list[1] = B;

    // // print_m(list[0]);
    // struct matrix *C = diag(list, 2);

    

    Quaternion *Q = NewQuarternion(1,2,3,4);
    quatDisplay(Q);



    
    return 0;
}