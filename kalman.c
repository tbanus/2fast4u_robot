
#define Nsta 7     // Two state values: pressure, temperature
#define Mobs 3     // Three measurements: baro pressure, baro temperature, LM35 temperature
#define timestep 0.002 //miliseconds
#define acc_timestep 0.000002 // khjk
#define enc_uncertainty 1.5696//19
#define acc_uncertainty 1.5696
#define gyro_uncertainty 1.5696
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#define timestep 0.002 //miliseconds
#define acc_timestep 0.000002 // khjk


typedef struct {

    int n;          /* number of state values */
    int m;          /* number of observables */

    double x[Nsta][1];    /* state vector */
    double z[Mobs][1];

    double P[Nsta][Nsta];  /* prediction error covariance */
    double Q[Nsta][Nsta];  /* process noise covariance */
    double R[Mobs][Mobs];  /* measurement error covariance */

    double G[Nsta][Mobs];  /* Kalman gain; a.k.a. K */

    double F[Nsta][Nsta];  /* Jacobian of process model */
    double H[Mobs][Nsta];  /* Jacobian of measurement model */

    double Ht[Nsta][Mobs]; /* transpose of measurement Jacobian */
    double Ft[Nsta][Nsta]; /* transpose of process Jacobian */
    double Gt[Mobs][Nsta];  /* Kalman gain; a.k.a. K */
    double Pp[Nsta][Nsta]; /* P, post-prediction, pre-update */

    double fx[Nsta][1];   /* output of user defined f() state-transition function */
    double hx[Mobs][1];   /* output of user defined h() measurement function */

    double I[Nsta][Nsta];  /* Unit Matrix */
    /* temporary storage */
    double tmp0[Nsta][Nsta];
    double tmp1[Nsta][Mobs];
    double tmp2[Mobs][Nsta];
    double tmp3[Mobs][Mobs];
    double tmp4[Mobs][Mobs];
    double tmp5[Mobs]; 
    double tmp6[Nsta][Nsta]; 
    double tmp7[Nsta][Nsta];
} ekf_t;        





double theta;
/* Cholesky-decomposition matrix-inversion code, adapated from
   http://jean-pierre.moreau.pagesperso-orange.fr/Cplus/choles_cpp.txt */
static void printf_matrix(double * matrix, int row_size,int col_size)
{
    int i, j, l;
    for(i=0; i < row_size ; ++i){
        for(j=0; j < col_size; ++j) 
        {    
            //printf("i=%d, j=%d, val= %f\t",i,j,matrix[i*col_size+j]);
            printf("%f\t",matrix[i*col_size+j]);
        }
        printf("\n");}
        printf("end\n");
}

static void eye(double * matrix, int size, double val)
{
    int i,j;
    for (i=0;i<size;++i)
    {
        for (j=0;j<size;++j)
        {

            if(i==j)
            {
                matrix[i*size+j]=val;
            }
            else
            {
                matrix[i*size+j]=0.0;
            }
        }
    }
}

static int choldc1(double * a, double * p, int n) {
    int i,j,k;
    double sum;

    for (i = 0; i < n; i++) {
        for (j = i; j < n; j++) {
            sum = a[i*n+j];
            for (k = i - 1; k >= 0; k--) {
                sum -= a[i*n+k] * a[j*n+k];
            }
            if (i == j) {
                if (sum <= 0) {
                    return 1; /* error */
                }
                p[i] = sqrt(sum);
            }
            else {
                a[j*n+i] = sum / p[i];
            }
        }
    }

    return 0; /* success */
}

static int choldcsl(double * A, double * a, double * p, int n) 
{
    int i,j,k; double sum;
    for (i = 0; i < n; i++) 
        for (j = 0; j < n; j++) 
            a[i*n+j] = A[i*n+j];
    if (choldc1(a, p, n)) return 1;
    for (i = 0; i < n; i++) {
        a[i*n+i] = 1 / p[i];
        for (j = i + 1; j < n; j++) {
            sum = 0;
            for (k = i; k < j; k++) {
                sum -= a[j*n+k] * a[k*n+i];
            }
            a[j*n+i] = sum / p[j];
        }
    }

    return 0; /* success */
}


static int cholsl(double * A, double * a, double * p, int n) 
{
    int i,j,k;
    if (choldcsl(A,a,p,n)) return 1;
    for (i = 0; i < n; i++) {
        for (j = i + 1; j < n; j++) {
            a[i*n+j] = 0.0;
        }
    }
    for (i = 0; i < n; i++) {
        a[i*n+i] *= a[i*n+i];
        for (k = i + 1; k < n; k++) {
            a[i*n+i] += a[k*n+i] * a[k*n+i];
        }
        for (j = i + 1; j < n; j++) {
            for (k = j; k < n; k++) {
                a[i*n+j] += a[k*n+i] * a[k*n+j];
            }
        }
    }
    for (i = 0; i < n; i++) {
        for (j = 0; j < i; j++) {
            a[i*n+j] = a[j*n+i];
        }
    }

    return 0; /* success */
}

static void zeros(double * a, int m, int n)
{
    int j;
    for (j=0; j<m*n; ++j)
        a[j] = 0;
}

#ifdef DEBUG
static void dump(double * a, int m, int n, const char * fmt)
{
    int i,j;

    char f[100];
    sprintf(f, "%s ", fmt);
    for(i=0; i<m; ++i) {
        for(j=0; j<n; ++j)
            printf(f, a[i*n+j]);
        printf("\n");
    }
}
#endif

/* C <- A * B */
static void mulmat(double * a, double * b, double * c, int arows, int acols, int bcols)
{
    int i, j,l;

    for(i=0; i<arows; ++i)
        for(j=0; j<bcols; ++j) {
            c[i*bcols+j] = 0;
            for(l=0; l<acols; ++l)
                c[i*bcols+j] += a[i*acols+l] * b[l*bcols+j];
        }
}

static void mulvec(double * a, double * x, double * y, int m, int n)
{
    int i, j;

    for(i=0; i<m; ++i) {
        y[i] = 0;
        for(j=0; j<n; ++j)
            y[i] += x[j] * a[i*n+j];
    }
}

static void transpose(double * a, double * at, int m, int n)
{
    int i,j;

    for(i=0; i<m; ++i)
        for(j=0; j<n; ++j) {
            at[j*m+i] = a[i*n+j];
        }
}

/* A <- A + B */
static void accum(double * a, double * b, int m, int n)
{        
    int i,j;

    for(i=0; i<m; ++i)
        for(j=0; j<n; ++j)
            a[i*n+j] += b[i*n+j];
}

/* C <- A + B */
static void add(double * a, double * b, double * c, int n)
{
    int j;

    for(j=0; j<n; ++j)
        c[j] = a[j] + b[j];
}


/* C <- A - B */
static void sub(double * a, double * b, double * c, int n)
{
    int j;

    for(j=0; j<n; ++j){
        c[j] = a[j] - b[j];
        //printf("Line is: %f\n",a[j]);
    }
}

static void negate(double * a, int m, int n)
{        
    int i, j;

    for(i=0; i<m; ++i)
        for(j=0; j<n; ++j)
            a[i*n+j] = -a[i*n+j];
}

static void mat_addeye(double * a, int n)
{
    int i;
    for (i=0; i<n; ++i)
        a[i*n+i] += 1;
}

static void copy_matrix(double * matrix, double * target, int row_size,int col_size)
{
    int i, j;

    for(i=0; i < row_size ; ++i)
    {
        for(j=0; j < col_size; ++j) 
        {   
            target[i*col_size+j]=matrix[i*col_size+j];
        }
    }
}


int kalman_init(ekf_t * ekf)
{


double H[Mobs][Nsta] = {
                        {0, 0 ,0, 1, 0, 0, 0}, 
                        {0, 0 ,0, 0, 1, 0, 0},
                        {0, 0 ,0, 0, 0, 1, 0} }; //Map measurement from  IMU 

copy_matrix(H, ekf->H,3,7);



eye(ekf->Q ,  Nsta, 0.01);

double x[Nsta][1] = {
                            {0},
                            {0},
                            {0},
                            {0},
                            {0},
                            {0},
                            {0},
                            };
copy_matrix(x, ekf->x,7,1);

theta=ekf->x[2][0];
double F[ Nsta][ Nsta] = {
            {1, 0, 0, timestep*cos(theta), 0, acc_timestep*cos(theta), 0},
            {0, 1, 0, timestep*sin(theta), 0, acc_timestep*sin(theta), 0},
            {0, 0, 1.0, 0, timestep, 0, acc_timestep},
            {0, 0, 0, 1.0, 0, timestep, 0},
            {0, 0, 0, 0, 1.0, 0, timestep},
            {0, 0, 0, 0, 0, 1.0, 0},
            {0, 0, 0, 0, 0, 0, 1.0}
            };
copy_matrix(F, ekf->F,7,7);

eye(ekf->P,7,5);


double R[Mobs][Mobs]=  {
                        {acc_uncertainty, 0, 0},
                        {0, acc_uncertainty, 0},
                        {0, 0, gyro_uncertainty}
                        }; 
                         /* measurement error covariance */
copy_matrix(R, ekf->R,3,3);

double I[Nsta][Nsta]; eye(I, Nsta, 1.0);
copy_matrix(I, ekf->I,7,7);







//State Extrapolation
mulmat(ekf->F, x,  ekf->fx,  Nsta,  Nsta, 1);

//Covariance Extrapolation
mulmat(ekf->F, ekf->P, ekf->tmp0, Nsta, Nsta, Nsta);

transpose(ekf->F, ekf->Ft, Nsta, Nsta);
mulmat(ekf->tmp0, ekf->Ft, ekf->Pp, Nsta, Nsta, Nsta);
accum(ekf->Pp, ekf->Q, Nsta, Nsta);



return 1 ;
}
//END OF INIT

int kalman_update(ekf_t* ekf){

double theta=ekf->x[2][0];

ekf->F[0][3]=timestep*cos(theta);
ekf->F[0][5]=acc_timestep*cos(theta);
ekf->F[1][3]=timestep*sin(theta);
ekf->F[1][5]=acc_timestep*sin(theta);





/* G_k = P_k H^T_k (H_k P_k H^T_k + R)^{-1} */
//CALC KALMAN GAIN
transpose(ekf->H,  ekf->Ht, Mobs, Nsta);
mulmat( ekf->Pp,  ekf->Ht,  ekf->tmp1, Nsta, Nsta, Mobs);
mulmat(ekf->H,  ekf->Pp,  ekf->tmp2, Mobs,  Nsta ,Nsta);
mulmat( ekf->tmp2,  ekf->Ht,  ekf->tmp3, Mobs,  Nsta , Mobs);
accum( ekf->tmp3,  ekf->R, Mobs,  Mobs);
if (cholsl( ekf->tmp3,  ekf->tmp4,  ekf->tmp5, Mobs)){printf("err"); return 1;}
mulmat( ekf->tmp1,  ekf->tmp4,  ekf->G, Nsta, Mobs,  Mobs);
mulmat(ekf->H,  ekf->fx, ekf->hx,  Mobs,  Nsta, 1);
sub(ekf->z,  ekf->hx,  ekf->tmp5, 3);
mulmat( ekf->G,  ekf->tmp5,  ekf->tmp2, Nsta, Mobs, 1);
add(  ekf->fx,  ekf->tmp2,  ekf->x, Nsta);
mulmat( ekf->G,  ekf->H,  ekf->tmp0, Nsta, Mobs,  Nsta);
sub(ekf->I , ekf->tmp0 , ekf->tmp6 ,49); 
mulmat(ekf->tmp6,ekf->Pp,ekf->tmp0,Nsta ,Nsta, Nsta);
transpose(ekf->tmp6,ekf->tmp7,Nsta ,Nsta);
mulmat(ekf->tmp0,ekf->tmp7, ekf->tmp6, Nsta,Nsta,Nsta); 
transpose(ekf->G,ekf->Gt,Nsta,Mobs);
mulmat(ekf->G,ekf->R,ekf->tmp1,Nsta,Mobs, Mobs);
mulmat(ekf->tmp1,ekf->Gt,ekf->tmp0,Nsta,Mobs, Nsta);
add(ekf->tmp0,ekf->tmp6,ekf->P,49);   

//State Extrapolation
mulmat(ekf->F, ekf->x,  ekf->fx,  Nsta,  Nsta, 1);

//Covariance Extrapolation
mulmat(ekf->F, ekf->P, ekf->tmp0, Nsta, Nsta, Nsta);
transpose(ekf->F, ekf->Ft, Nsta, Nsta);
mulmat(ekf->tmp0, ekf->Ft, ekf->Pp, Nsta, Nsta, Nsta);
accum(ekf->Pp, ekf->Q, Nsta, Nsta);




 // #TODO: Clear temps
 // #TODO: use *ekf->G 
 // #TODO: use QR Update




return 1;

}
