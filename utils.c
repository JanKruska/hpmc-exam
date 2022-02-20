#include <stdio.h>
#include <errno.h>  // for errno
#include <limits.h> // for INT_MAX, INT_MIN
#include <stdlib.h> // for strtol
#include <math.h>
#include "omp.h"

#define A(i, j) *(ap + (j)*lda + (i)) // map A( i,j )    to array ap    in column-major order

void RandomMatrix(int m, int n, double *ap, int lda)
/*
   RandomMatrix overwrite A with random values.
*/
{
  int i, j;

#pragma omp parallel for private(i, j)
  for (j = 0; j < n; j++)
  {
    for (i = 0; i < m; i++)
    {
      A(i, j) = drand48();
    }
  }
}

void ZeroMatrix(int m, int n, double *ap, int lda)
/*
   RandomMatrix overwrite A with random values.
*/
{
  int i, j;

#pragma omp parallel for private(i, j)
  for (j = 0; j < n; j++)
  {
    for (i = 0; i < m; i++)
    {
      A(i, j) = 0.0;
    }
  }
}

void RandomMatrixSymmetric(int n, double *ap, int lda)
/*
   RandomMatrix Symmetric Positive definite.
   #TODO: Only generate upper(or lower) triangular matrix since BLAS ignores the rest
*/
{
  RandomMatrix(n, n, ap, lda);
  int i, j;
  double a, b;

#pragma omp parallel for private(i, j)
  for (j = 0; j < n; j++)
  {
    for (i = j + 1; i < n; i++)
    {
      a = A(i, j);
      b = A(j, i);
      A(i, j) = A(j, i) = 0.5 * (a + b);
    }
  }
#pragma omp parallel for private(i, j)
  for (i = 0; i < n; i++)
  {
    A(i, i) = A(i, i) + n;
  }
  // #pragma omp parallel for private(i,j)
  // for ( j=0; j<n; j++ ) {
  //   for ( i=j+1; i<n; i++ ) {
  //     A( j,i ) = A( i,j );
  //   }
  // }
}

void printMatrix(int m, int n, double *ap, int lda)
{
  int i, j;
  for (i = 0; i < n; i++)
  {
    for (j = 0; j < m; j++)
    {
      printf("%8.2f ", A(i, j));
    }
    printf("\n");
  }
}

int str2int(char *string)
{
  char *p;
  int num;

  errno = 0;
  long conv = strtol(string, &p, 10);

  // Check for errors: e.g., the string does not represent an integer
  // or the integer is larger than int
  if (errno != 0 || *p != '\0' || conv > INT_MAX || conv < INT_MIN)
  {
    // Put here the handling of the error, like exiting the program with
    // an error message
  }
  else
  {
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
  int i, j;
  if (uplo == 'u')
  {
#pragma omp parallel for private(i, j)
    for (j = 0; j < n; j++)
    {
      for (i = j + 1; i < n; i++)
      {
        A(i, j) = A(j, i);
      }
    }
  }
  else if (uplo == 'l')
  {
#pragma omp parallel for private(i, j)
    for (j = 0; j < n; j++)
    {
      for (i = j + 1; i < n; i++)
      {
        A(j, i) = A(i, j);
      }
    }
  }
}

int issymetric(double *ap, int n, int m, int lda)
{
  int i, j;

  // #pragma omp parallel for private(i,j)
  for (j = 0; j < n; j++)
  {
    for (i = j + 1; i < n; i++)
    {
      if (A(i, j) != A(j, i))
      {
        return 0;
      }
    }
  }
  return 1;
}

double arrayMin(double *x, int size)
{
  int i;
  double min = x[0];
  for (i = 1; i < size; i++)
  {
    min = (x[i] < min ? x[i] : min);
  }
  return min;
}

double arrayMax(double *x, int size)
{
  int i;
  double max = x[0];
  for (i = 1; i < size; i++)
  {
    max = (x[i] > max ? x[i] : max);
  }
  return max;
}

double arraySum(double *x, int size)
{
  int i;
  double sum = 0.0;
  for (i = 0; i < size; i++)
  {
    sum += x[i];
  }
  return sum;
}

double arrayMean(double *x, int size)
{
  return (arraySum(x, size) / size);
}

double arrayVariance(double *x, int size)
{
  int i;
  double mean = arrayMean(x, size), sum = 0.0;
  for (i = 0; i < size; i++)
  {
    sum += (x[i] - mean) * (x[i] - mean);
  }
  return sum / size;
}

double arrayStd(double *x, int size)
{
  return (sqrt(arrayVariance(x, size)));
}
