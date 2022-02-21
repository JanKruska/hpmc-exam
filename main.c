#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <getopt.h>
#include "main.h"
#include "utils.h"
#include "omp.h"

struct result
{
  // double *matrix;
  double time;
  double gflops;
};

/* Flag set by ‘--parsable. Exports output as csv formatted text to stdout*/
static int parsable_flag;
/* Flag set by ‘--no-header. */
static int no_header_flag;

void usage()
{
  printf("Usage: ./main [options] M \n"
         "or:    ./main [options] l u s \n"
         "M: matrix dimension m of problem; l,u,s: lower, upper and step for matrix dimension m");
  exit(EXIT_FAILURE);
}

struct result dcase1(int m, double *C)
{
  int n, k,
      ldA, ldB, ldY, ldC,
      size, first, last, inc;

  double
      d_one = 1.0,
      d_zero = 0.0,
      d_half = 0.5,
      dtime, tic, gflops;

  double
      *A, *B, *Y;

  first = 120;
  last = 300;
  inc = 20;
  /* Adjust first and last so that they are multiples of inc */
  last = (last / inc) * inc;
  first = (first / inc) * inc;
  first = (first == 0 ? inc : first);

  ldC = m;
  
  ZeroMatrix(m, m, C, ldC);
  // dtime = omp_get_wtime();
  dtime = gflops = 0.0;
  for (n = last; n >= first; n -= inc)
  {
    k = ldA = ldB = ldY = n;

    /* Gflops performed */
    /*FLOPS taken from http://www.netlib.org/lapack/lawnspdf/lawn41.pdf */
    gflops += 1e-09 * (2.0 * m * m * n + 2.0 * m * n * n + n);

    /* Allocate space for the matrices. */
    /* A n*n, B n*m, Y n*m, C m*m*/
    A = (double *)malloc(ldA * ldB * sizeof(double));
    B = (double *)malloc(ldB * m * sizeof(double));
    Y = (double *)malloc(ldA * m * sizeof(double));

    /* Generate random matrix A */
    RandomMatrixSymmetric(n, A, ldA);

    /* Generate random matrix B */
    RandomMatrix(n, m, B, ldB);

    tic = omp_get_wtime();

    dsymm_("l", "u",
           &n, &m,
           &d_half, A, &ldA,
           B, &ldB,
           &d_zero, Y, &ldY);

    dsyr2k_("u", "t",
            &m, &n,
            &d_one, B, &ldB,
            Y, &ldY,
            &d_one, C, &ldC);

    /* stop clock */
    dtime += omp_get_wtime() - tic;

    /* Free the buffers */
    free(A);
    free(B);
    free(Y);
  }
  resymmetrize(m, m, C, ldC, 'u');

  struct result r;
  r.time = dtime;
  r.gflops = gflops;
  // r.matrix = C;
  return r;
}

int main(int argc, char *argv[])
{
  int m, min_m, max_m, step_m,
      irep, nrepeats,
      c;

  double
      *C, *times, *array_gflops;

  nrepeats = 1;
  while (1)
  {
    static struct option long_options[] =
        {
            /* These options set a flag. */
            {"parsable", no_argument, &parsable_flag, 1},
            {"no-header", no_argument, &no_header_flag, 1},
            /* These options don’t set a flag.
               We distinguish them by their indices. */
            {"repetitions", required_argument, 0, 'r'},
            {0, 0, 0, 0}};
    /* getopt_long stores the option index here. */
    int option_index = 0;

    c = getopt_long(argc, argv, ":r:",
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
      printf("option %s", long_options[option_index].name);
      if (optarg)
        printf(" with arg %s", optarg);
      printf("\n");
      break;

    case 'r':
      nrepeats = str2int(optarg);
      // printf ("option -d with value `%s'\n", optarg);
      break;

    case '?':
      /* getopt_long already printed an error message. */
      break;

    default:
      abort();
    }
  }

  if (argc - optind == 1)
  {
    min_m = max_m = str2int(argv[optind]);
    step_m = 1.0;
  }
  else if (argc - optind == 3)
  {
    min_m = str2int(argv[optind]);
    max_m = str2int(argv[optind + 1]);
    step_m = str2int(argv[optind + 2]);
  }
  else
  {
    usage();
  }

  if (!parsable_flag)
    printf("m=%d to %d with increment %d for %d repetitions\n", min_m, max_m, step_m, nrepeats);

  if (!no_header_flag)
    if (parsable_flag)
      printf("m,time,time_std,GFLOPS,GFLOPS_mean,GFLOPS_std,GFLOPS/core\n");
    else
      printf("   m     time      time_std   \t GFLOPS GFLOPS_mean GFLOPS_std GFLOPS/core\n");

  times = malloc(nrepeats * sizeof *times);
  array_gflops = malloc(nrepeats * sizeof *array_gflops);
  struct result r;
  for (m = min_m; m <= max_m; m += step_m)
  {
    C = (double *)malloc(m * m * sizeof(double));
    for (irep = 0; irep < nrepeats; irep++)
    {
      r = dcase1(m,C);
      /* record the time */
      times[irep] = r.time;
      array_gflops[irep] = r.gflops / r.time;
    }
    free(C);
    if (parsable_flag)
    {
      printf("%d,%e,%e,%e,%e,%e,%e\n",
             m,
             arrayMin(times, nrepeats),
             arrayStd(times, nrepeats),
             arrayMax(array_gflops, nrepeats),
             arrayMean(array_gflops, nrepeats),
             arrayStd(array_gflops, nrepeats),
             arrayMax(array_gflops, nrepeats) / omp_get_max_threads());
    }
    else
    {
      printf(" %5d %8.4le %8.4le %10.4f %10.4f %10.4f %10.4f\n",
             m,
             arrayMin(times, nrepeats),
             arrayStd(times, nrepeats),
             arrayMax(array_gflops, nrepeats),
             arrayMean(array_gflops, nrepeats),
             arrayStd(array_gflops, nrepeats),
             arrayMax(array_gflops, nrepeats) / omp_get_max_threads());
    }
    // We flush the output buffer because otherwise
    // it may throw the timings of a next
    // experiment.
    fflush(stdout);
  }
  free(times);
  free(array_gflops);
  exit(0);
}
