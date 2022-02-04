#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <errno.h>   // for errno
#include <limits.h>  // for INT_MAX, INT_MIN
#include <stdlib.h>  // for strtol
#include "omp.h"

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

#define A( i,j ) *( ap + (j)*lda + (i) )          // map A( i,j )    to array ap    in column-major order

void RandomMatrix( int m, int n, double *ap, int lda )
/* 
   RandomMatrix overwrite A with random values.
*/
{
  int  i, j;
  
  #pragma omp parallel for private(i,j)
  for ( j=0; j<n; j++ ) {
    for ( i=0; i<m; i++ ) {
      A( i,j ) = drand48();
    }
  }
}

void ZeroMatrix( int m, int n, double *ap, int lda )
/* 
   RandomMatrix overwrite A with random values.
*/
{
  int  i, j;
  
  #pragma omp parallel for private(i,j)
  for ( j=0; j<n; j++ ) {
    for ( i=0; i<m; i++ ) {
      A( i,j ) = 0.0;
    }
  }
}

void RandomMatrixSymmetric(int n, double *ap, int lda )
/* 
   RandomMatrix Symmetric Positive definite.
   #TODO: Only generate upper(or lower) triangular matrix since BLAS ignores the rest
*/
{
  RandomMatrix(n,n,ap,lda);
  int  i, j;
  double a,b;
  
  #pragma omp parallel for private(i,j)
  for ( j=0; j<n; j++ ) {
    for ( i=j+1; i<n; i++ ) {
      a = A( i,j ) ;
      b = A( j,i ) ;
      A( i,j ) = A( j,i ) = 0.5*(a+b);
    }
  }
  #pragma omp parallel for private(i,j)
  for ( i=0; i<n; i++ ) {
      A( i,i ) = A( i,i ) + n;
  }
  // #pragma omp parallel for private(i,j)
  // for ( j=0; j<n; j++ ) {
  //   for ( i=j+1; i<n; i++ ) {
  //     A( j,i ) = A( i,j );
  //   }
  // }
}

void printMatrix(int m, int n, double *ap, int lda){
  int i,j;
  for(i = 0; i < n; i++) {
        for(j = 0; j < m; j++) {
            printf("%8.2f ", A( i,j ));
        }
        printf("\n");
    } 
}

void usage(){
    printf("Usage: ./main M \n M - matrix dimension m of problem");
    exit(EXIT_FAILURE);
}

int str2int(char * string){
    char *p;
    int num;

    errno = 0;
    long conv = strtol(string, &p, 10);

    // Check for errors: e.g., the string does not represent an integer
    // or the integer is larger than int
    if (errno != 0 || *p != '\0' || conv > INT_MAX || conv < INT_MIN) {
        // Put here the handling of the error, like exiting the program with
        // an error message
    } else {
        num = conv;
        return num;
    }
}

void resymmetrize(int m, int n, double *ap, int lda, char uplo)
/* 
   RandomMatrix Symmetric Positive definite.
   #TODO: Only generate upper(or lower) triangular matrix since BLAS ignores the rest
*/
{
  int  i, j;
  if(uplo=='u'){
    #pragma omp parallel for private(i,j)
    for ( j=0; j<n; j++ ) {
      for ( i=j+1; i<n; i++ ) {
        A( i,j ) = A( j,i );
      }
    }
  }
  else if(uplo=='l'){
    #pragma omp parallel for private(i,j)
    for ( j=0; j<n; j++ ) {
      for ( i=j+1; i<n; i++ ) {
        A( j,i ) = A( i,j );
      }
    }
  }
}

int issymetric(double *ap, int n, int m, int lda){
  int  i, j;
  
  // #pragma omp parallel for private(i,j)
  for ( j=0; j<n; j++ ) {
    for ( i=j+1; i<n; i++ ) {
      if(A(i,j) != A(j,i)){
        return 0;
      }
    }
  }
  return 1;
}

int main(int argc, char *argv[])
{
  int
    x, y,
    m, n, k,
    ldA, ldB, ldY, ldC,
    size, first, last, inc,
    i, irep,
    nrepeats;

  double
    d_one = 1.0,
    d_zero = 0.0,
    d_half = 0.5,
    dtime, dtime_best, 
    diff, maxdiff = 0.0, gflops;

  double
    *A, *B, *Y, *C;

  /* Print the number of threads available */
//   printf( "%% Number of threads = %d\n\n", omp_get_max_threads() );
  /* Every time trial is repeated "repeat" times and the fastest run in recorded */
//   printf( "%% size m of matrix:" );
//   scanf( "%d", &m );
//   printf( "%% %d\n", m );
if(argc<=1 || argc > 2){
    usage();
}else{
    m = str2int(argv[1]);
}
    nrepeats = 1;
  /* Timing trials for matrix sizes m=n=k=first to last in increments
     of inc will be performed.  (Actually, we are going to go from
     largest to smallest since this seems to give more reliable 
     timings.  */
//   printf( "%% enter first, last, inc:" );
//   scanf( "%d%d%d", &first, &last, &inc );
//  int n = 100
first = 120;
last = 300;
inc = 20;
  /* Adjust first and last so that they are multiples of inc */
  last = ( last / inc ) * inc;
  first = ( first / inc ) * inc;
  first = ( first == 0 ? inc : first );
  
  printf( "n = %d to %d with increment %d, m=%d\n", first, last, inc, m );

  printf( "data = [\n" );
  printf( "%%  n     time       GFLOPS  GFLOPS/core\n" );
  
  C = ( double * ) malloc( ldC * m * sizeof( double ) );
  ldC = m;
  /* Generate random matrix C */
  ZeroMatrix( m, m, C, ldC );

  for ( size=last; size>= first; size-=inc ){
    /* we will only time cases where all three matrices are square */
    n = k = size;
    ldA = ldB = ldY = size;
    
    /* Gflops performed */
    gflops = 2.0 * m * n * k * 1e-09;

    /* Allocate space for the matrices. */
    /* A n*n, B n*m, Y n*m, C m*m*/
    A = ( double * ) malloc( ldA * ldB * sizeof( double ) );
    B = ( double * ) malloc( ldB * m * sizeof( double ) );
    Y = ( double * ) malloc( ldA * m * sizeof( double ) );

    /* Generate random matrix A */
    RandomMatrixSymmetric(n, A, ldA );

    /* Generate random matrix B */
    RandomMatrix( n, m, B, ldB );

    /* Time dgemm (double precision general matrix-matrix
       multiplicationn */
    for ( irep=0; irep<nrepeats; irep++ ){
    
      /* start clock */
      dtime = omp_get_wtime();

      dsymm_( "l", "u",
        &n, &m,
        &d_half, A, &ldA,
                B, &ldB,
        &d_zero, Y, &ldY );

      dsyr2k_( "u", "t",
        &m, &n,
        &d_one, B, &ldB,
                Y, &ldY,
        &d_one, C, &ldC );

      /* stop clock */
      dtime = omp_get_wtime() - dtime;

      /* record the best time so far */
      if ( irep == 0 ) 
        dtime_best = dtime;
      else
        dtime_best = ( dtime < dtime_best ? dtime : dtime_best );
    }
  
//     printf( " %5d %8.4le %8.4le %8.4le\n", n, dtime_best, gflops/dtime_best, gflops/dtime_best/omp_get_max_threads() );
    printf( " %5d %8.4le %8.4f %8.4f\n", n, dtime_best, gflops/dtime_best, gflops/dtime_best/omp_get_max_threads() );
    
    // printMatrix(m,n,C,ldY);
    // printMatrix(n,n,A,ldA);
    // printMatrix(m,n,C,ldY);
    
    
    fflush( stdout );  // We flush the output buffer because otherwise
		       // it may throw the timings of a next
		       // experiment.
    

    /* Free the buffers */
    free( A );
    free( B );
    free( Y );
  }
  resymmetrize(m,m,C,ldC,'u');
  printMatrix(m,m,C,ldC);
  free( C );
  printf( "];\n\n" );
  
  exit( 0 );
}
