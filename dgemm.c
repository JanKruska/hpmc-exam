#include <stdio.h>
#include <stdlib.h>
#include "omp.h"
#include "mkl.h"
#include "utils.h"
 
int main(void)
  {
    MKLVersion Version;
 
    mkl_get_version(&Version);
 
 
    printf("Major version:           %d\n",Version.MajorVersion);
    printf("Minor version:           %d\n",Version.MinorVersion);
    printf("Update version:          %d\n",Version.UpdateVersion);
    printf("Product status:          %s\n",Version.ProductStatus);
    printf("Build:                   %s\n",Version.Build);
    printf("Platform:                %s\n",Version.Platform);
    printf("Processor optimization:  %s\n",Version.Processor);
    printf("================================================================\n");
    printf("\n");


    int
    x, y,
    m, n, k,
    ldA, ldB, ldC,
    i, irep,
    nrepeats;

  double
    d_one = 1.0,
    d_zero = 0.0,
    dtime, dtime_best, dtime_std, tic,
    diff, maxdiff = 0.0, gflops;

  double
    *A, *B, *C;
    nrepeats = 10;

  double *times = malloc(nrepeats * sizeof *times);    
double *array_gflops = malloc(nrepeats * sizeof *times);  
n = m = k = 1000;
printf( "   m     time      time_std   \t GFLOPS GFLOPS_mean GFLOPS_std GFLOPS/core\n" );
for ( irep=0; irep<nrepeats; irep++ ){
  /* start clock */
  // dtime = omp_get_wtime();
  dtime = gflops = 0.0;
   ldA = ldB = ldC = n;
    
    /* Gflops performed */
    /*FLOPS taken from http://www.netlib.org/lapack/lawnspdf/lawn41.pdf */
    gflops += 1e-09 * (2.0*m*k*n);

    A = ( double * ) malloc( ldA * ldB * sizeof( double ) );
    B = ( double * ) malloc( ldB * m * sizeof( double ) );
    C = ( double * ) malloc( ldB * m * sizeof( double ) );

    /* Generate random matrix A */
    RandomMatrix(n,m, A, ldA );

    /* Generate random matrix B */
    RandomMatrix( m,k, B, ldB );
    
    tic = omp_get_wtime();

    dgemm_( "n", "n",
	      &n, &m, &n,
	      &d_one, A, &ldA,
	              B, &ldB,
	      &d_zero, C, &ldC );

    /* stop clock */
    dtime += omp_get_wtime() - tic;
    
    /* Free the buffers */
    free( A );
    free( B );
    free( C );
    
    /* record the time */
    times[irep] = dtime;
    array_gflops[irep] = gflops/dtime;
  
  // printMatrix(m,m,C,ldC);
  }
  dtime_best = arrayMin(times,nrepeats);
  dtime_std = arrayStd(times,nrepeats);
  if(0){
    printf( "%d,%e,%e,%e,%e,%e,%e\n", 
    m, 
    arrayMin(times,nrepeats), 
    arrayStd(times,nrepeats), 
    arrayMax(array_gflops,nrepeats),
    arrayMean(array_gflops,nrepeats),
    arrayStd(array_gflops,nrepeats), 
    arrayMax(array_gflops,nrepeats)/omp_get_max_threads() );
  }else{
    printf( " %5d %8.4le %8.4le %10.4f %10.4f %10.4f %10.4f\n", 
    m, 
    arrayMin(times,nrepeats), 
    arrayStd(times,nrepeats), 
    arrayMax(array_gflops,nrepeats),
    arrayMean(array_gflops,nrepeats),
    arrayStd(array_gflops,nrepeats), 
    arrayMax(array_gflops,nrepeats)/omp_get_max_threads() );
  }

    return 0;
  }
