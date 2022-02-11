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