#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <errno.h>   // for errno
#include <limits.h>  // for INT_MAX, INT_MIN
#include <stdlib.h>  // for strtol
#include "omp.h"
#include <getopt.h>
#include "main.h"

/* Flag set by ‘--parsable. */
static int parsable_flag;
/* Flag set by ‘--no-header. */
static int no_header_flag;

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
    min_m,max_m,step_m,
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
// if(argc == 2){
//     m = str2int(argv[1]);
// }else{
//     usage();
// }
int aflag = 0;
  int bflag = 0;
  char *cvalue = NULL;
  int index;
  int c;
  while (1)
    {
      static struct option long_options[] =
        {
          /* These options set a flag. */
          {"parsable", no_argument,       &parsable_flag, 1},
          {"no-header", no_argument,       &no_header_flag, 1},
          // {"brief",   no_argument,       &verbose_flag, 0},
          /* These options don’t set a flag.
             We distinguish them by their indices. */
          {"add",     no_argument,       0, 'a'},
          {"append",  no_argument,       0, 'b'},
          {"delete",  required_argument, 0, 'd'},
          {"create",  required_argument, 0, 'c'},
          {"file",    required_argument, 0, 'f'},
          {0, 0, 0, 0}
        };
      /* getopt_long stores the option index here. */
      int option_index = 0;

      c = getopt_long (argc, argv, "abc:d:f:",
                       long_options, &option_index);

      /* Detect the end of the options. */
      if (c == -1)
        break;

      switch (c)
        {
        case 0:
          /* If this option set a flag, do nothing else now. */
          if (long_options[option_index].flag != 0)
            break;
          printf ("option %s", long_options[option_index].name);
          if (optarg)
            printf (" with arg %s", optarg);
          printf ("\n");
          break;

        case 'a':
          puts ("option -a\n");
          break;

        case 'b':
          puts ("option -b\n");
          break;

        case 'c':
          printf ("option -c with value `%s'\n", optarg);
          break;

        case 'd':
          printf ("option -d with value `%s'\n", optarg);
          break;

        case 'f':
          printf ("option -f with value `%s'\n", optarg);
          break;

        case '?':
          /* getopt_long already printed an error message. */
          break;

        default:
          abort ();
        }
    }


  // printf ("aflag = %d, bflag = %d, cvalue = %s\n",
  //         aflag, bflag, cvalue);
  // printf("%d %d \n",optind,argc);
  if(argc-optind==1){
    min_m = max_m = str2int(argv[optind]);
    step_m = 1.0;
  }else if(argc-optind==3){
    min_m = str2int(argv[optind]);
    max_m = str2int(argv[optind+1]);
    step_m = str2int(argv[optind+2]);
  }else
  {
    usage();
  }
  // for (index = optind; index < argc; index++)
    

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
  
  if(!parsable_flag)
    printf( "n = %d to %d with increment %d, m=%d to %d with increment %d\n", first, last, inc, min_m, max_m, step_m );

  if(!no_header_flag)
    printf( "   m     time       GFLOPS  GFLOPS/core\n" );
  
for(m=min_m; m<=max_m; m+=step_m){
ldC = m;
C = ( double * ) malloc( ldC * m * sizeof( double ) );
  
for ( irep=0; irep<nrepeats; irep++ ){
  ZeroMatrix( m, m, C, ldC );
  /* start clock */
  dtime = omp_get_wtime();
  gflops = 0.0;
  for ( n=last; n>= first; n-=inc ){
    k = ldA = ldB = ldY = n;
    
    /* Gflops performed */
    /*FLOPS taken from http://www.netlib.org/lapack/lawnspdf/lawn41.pdf */
    gflops += 1e-09 * (2.0*m*m*n + 2.0*m*n*n+n);

    /* Allocate space for the matrices. */
    /* A n*n, B n*m, Y n*m, C m*m*/
    A = ( double * ) malloc( ldA * ldB * sizeof( double ) );
    B = ( double * ) malloc( ldB * m * sizeof( double ) );
    Y = ( double * ) malloc( ldA * m * sizeof( double ) );

    /* Generate random matrix A */
    RandomMatrixSymmetric(n, A, ldA );

    /* Generate random matrix B */
    RandomMatrix( n, m, B, ldB );
    
    

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

    // printMatrix(m,n,C,ldY);
    // printMatrix(n,n,A,ldA);
    // printMatrix(m,n,C,ldY);
    
    }
    /* Free the buffers */
    free( A );
    free( B );
    free( Y );
     // We flush the output buffer because otherwise
            // it may throw the timings of a next
            // experiment.
    resymmetrize(m,m,C,ldC,'u');
    
    /* stop clock */
    dtime = omp_get_wtime() - dtime;

    /* record the best time so far */
    if ( irep == 0 ) 
      dtime_best = dtime;
    else
      dtime_best = ( dtime < dtime_best ? dtime : dtime_best );
  
  // printMatrix(m,m,C,ldC);
  }
//     printf( " %5d %8.4le %8.4le %8.4le\n", n, dtime_best, gflops/dtime_best, gflops/dtime_best/omp_get_max_threads() );
  printf( " %5d %8.4le %8.4f %8.4f\n", m, dtime_best, gflops/dtime_best, gflops/dtime_best/omp_get_max_threads() );
  fflush( stdout ); 
  free( C );
}
  exit( 0 );
}
