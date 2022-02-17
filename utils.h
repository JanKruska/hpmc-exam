void RandomMatrix( int m, int n, double *ap, int lda );
void ZeroMatrix( int m, int n, double *ap, int lda );
void RandomMatrixSymmetric(int n, double *ap, int lda );
void printMatrix(int m, int n, double *ap, int lda);
int str2int(char * string);
void resymmetrize(int m, int n, double *ap, int lda, char uplo);
int issymetric(double *ap, int n, int m, int lda);
double arrayMin(double *x, int size);
double arraySum(double *x, int size);
double arrayMean(double *x, int size);
double arrayStd(double *x, int size);
double arrayVariance(double *x, int size);
double arrayMax(double *x, int size);

/* Prototype for BLAS matrix-matrix multiplication routine (which we will 
   use for the reference implementation */
void dgemm_( char *, char *,                 // transA, transB
	     int *, int *, int *,            // m, n, k
	     double *, double *, int *,      // alpha, A, ldA
	               double *, int *,      //        B, ldB
	     double *, double *, int * );    // beta,  C, ldC

void dsymm_( char *, char *,                 // SIDE, UPLO
	     int *, int *,            // m, n
	     double *, double *, int *,      // alpha, A, ldA
	               double *, int *,      //        B, ldB
	     double *, double *, int * );    // beta,  C, ldC

void dsyr2k_( char *, char *,                 // UPLO, trans
	     int *, int *,            // m, n
	     double *, double *, int *,      // alpha, A, ldA
	               double *, int *,      //        B, ldB
	     double *, double *, int * );    // beta,  C, ldC

void dgemv_( char *,                     // trans,
	     int *, int *,                   // m, k
	     double *, double *, int *,      // alpha, A, ldA
	               double *, int *,      //        X, incx
	     double *, double *, int * );    // beta,  Y, incy