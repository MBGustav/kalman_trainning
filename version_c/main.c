#include <stdio.h>

#include "KFilter.h"
#include "structs.h"
#include "linear_algebra.h"


#define EPSILON (10.0f)
#define isEqual(x,y)     (fabs(fabs(x) - fabs(y)) <= EPSILON)

#ifndef FILEINPUT
#define FILEINPUT ("input.txt")
#endif //FILEINPUT



#ifndef NSAMPLES
#define NSAMPLES 100000
#endif

struct vector *GetRadar(double dT){ 
    static bool FirstRead = true;
    static int counter = 0;
    
    static struct vector** v_data;

    //create an array of vectors
    FILE *f; 
    // check_memory(DATA(z));
    if(FirstRead) {
        if((f = fopen(FILEINPUT, "r")) ==NULL){
            printf("Erro leitura Arquivo"); 
            exit(1);
        }
        v_data = (struct vector**)malloc(NSAMPLES * sizeof(struct vector*));
        //Leitura dado FILEINPUT
        float lect;
        for(int i = 0; i < NSAMPLES; i++){
            fscanf(f, "%f", &lect);
            fflush(stdin);
            v_data[i] = vector_new(1);
            VECTOR_IDX_INTO(v_data[i], 0) = lect;
        }
        FirstRead = false;
        fclose(f);
    //first iteration
    }

    //As long we filed vector, we should iterate over
    if(counter < NSAMPLES) 
        return v_data[counter++];
    //end iterations
    printf("Reached at the end of list - GetRadar(%d)\n", counter);
    exit(1);
}

bool
checkArrays(double *X, double *Y){
    FILE *Xf,*Zf; 
    if((Xf = fopen("XSaved.mtx", "r")) ==NULL)exit(1);
    if((Zf = fopen("ZSaved.mtx", "r"))==NULL)exit(1);
    bool equal;
    float alt, pos, vel;
    for(int i = 0; i < NSAMPLES; i++){
       fscanf(Xf, "%f;%f;%f",&pos, &vel, &alt);
       equal = isEqual(pos, X[3*(i)]);
       equal = equal && isEqual(vel, X[3*(i)+1]);
       equal = equal && isEqual(alt, X[3*(i)+2]);
        if(!equal){
            printf("Erro em X[%d]:\n", i);
            printf("\t   %5.5f => %5.5f\n", X[3*i]  , pos);
            printf("\t   %5.5f => %5.5f\n", X[3*(i)+1], vel);
            printf("\t %5.5f => %5.5f\n", X[3*(i)+2], alt);
            return false;
            
        }


    }
    printf("Vetores Com valores Próximos");
    return true;

}

int 
main(void){


//Set Initial Parameters:
    int m = 3, n = 1;
    double kappa = 0;     

    struct ss_d *_ss = (struct ss_d*) malloc(sizeof(struct ss_d));
    struct wm *_wm = (struct wm*) malloc(sizeof(struct wm));

    struct vector *x = vector_new(m); 
    VECTOR_IDX_INTO(x, 0) = 0.0;
    VECTOR_IDX_INTO(x, 1) = 90.0;
    VECTOR_IDX_INTO(x, 2) = 1100.0;
    _ss->x = x;

    struct matrix* R = matrix_new(n, n);
    MATRIX_IDX_INTO(R, 0,0) = 100.0;
    _wm->R = R;

    struct matrix* Q = matrix_identity(3); 
    matrix_scal_multiply(Q, Q, 0.01);
    _wm->Q = Q;

    struct matrix* P = matrix_identity(m);
    matrix_scal_multiply(P,P, 100.0);

    _wm->P = P;

    double *XSaved = (double*) malloc(NSAMPLES * sizeof(double)*3);
    double *ZSaved = (double*) malloc(NSAMPLES * sizeof(double));

    KFilter *KF = KF_Init(_ss, _wm, dt, kappa);

    struct matrix* A = matrix_new(3,3); 
    
    MATRIX_IDX_INTO(A, 0, 0) = 2;MATRIX_IDX_INTO(A, 0, 1) = 3;MATRIX_IDX_INTO(A, 0, 2) = 1;
    MATRIX_IDX_INTO(A, 1, 0) = 3;MATRIX_IDX_INTO(A, 1, 1) = 3;MATRIX_IDX_INTO(A, 1, 2) = 1;
    MATRIX_IDX_INTO(A, 2, 0) = 2;MATRIX_IDX_INTO(A, 2, 1) = 4;MATRIX_IDX_INTO(A, 2, 2) = 1;

    // B = matrix_inverse(A);
    // matrix_print(A);
    // matrix_print(B);
    // struct matrix*out = cholesky(A, 1.0);
    // matrix_print(out);
 
    struct vector* r = GetRadar(dt);
    
    printf(" ");

    double pos, vel, alt;
    for(int  itr =0; itr < NSAMPLES-1; itr++){
        
        r = GetRadar(dt);

        // printf("r[%d] = \n",itr);vector_print(r);
        KF_Update(KF, r);
        pos = VECTOR_IDX_INTO(x, 0);
        vel = VECTOR_IDX_INTO(x, 1);
        alt = VECTOR_IDX_INTO(x, 2);
        XSaved[3*(itr)+0] = pos; 
        XSaved[3*(itr)+1] = vel; 
        XSaved[3*(itr)+2] = alt; 
        // ZSaved[i] = VECTOR_IDX_INTO(r,0); 
        
    }
    #ifdef NDEBUG
    if(checkArrays(XSaved, ZSaved))
        printf("Execucao dentro da margem de Erro\n");
    else
        printf("Execucao com saida Divergente\n");

    
    printf("\nExecucao encerrada\n");

    #endif // NDEBUG


    return 0;
}
