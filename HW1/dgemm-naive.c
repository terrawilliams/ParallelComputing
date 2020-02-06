/* 
    Please include compiler name below (you may also include any other modules you would like to be loaded)

COMPILER= gnu

    Please include All compiler flags and libraries as you want them run. You can simply copy this over from the Makefile's first few lines
 
CC = cc
OPT = -O3
CFLAGS = -Wall -std=gnu99 $(OPT)
MKLROOT = /opt/intel/composer_xe_2013.1.117/mkl
LDLIBS = -lrt -Wl,--start-group $(MKLROOT)/lib/intel64/libmkl_intel_lp64.a $(MKLROOT)/lib/intel64/libmkl_sequential.a $(MKLROOT)/lib/intel64/libmkl_core.a -Wl,--end-group -lpthread -lm

*/

#include <stdio.h>
#include <smmintrin.h>
#include <immintrin.h>

const char* dgemm_desc = "Naive, three-loop dgemm.";

/* This routine performs a dgemm operation
 *  C := C + A * B
 * where A, B, and C are lda-by-lda matrices stored in column-major format.
 * On exit, A and B maintain their input values. */    
void square_dgemm (int n, double* A, double* B, double* C)
{
  /* For each row i of A */
  for (int j = 0; j < n; ++j)
    /* For each column j of B */
    for (int k = 0; k < n; k++)
    {
      for( int i = 0; i < n - 3; i += 4)
      {
          __m256d m1 = _mm256_load_pd(A + i + k * n);
          __m256d m2 = _mm256_set1_pd(*(B + k + j * n));
          __m256d m0 = _mm256_load_pd(C + i + j * n);
          m0 = _mm256_fmadd_pd(m1, m2, m0);
          _mm256_store_pd(C + i + j * n, m0);
      }
      if(n % 4 != 0)
      {
          if((n - 1) % 4 == 0)
          {
              C[(n - 1) + j * n] += A[(n - 1) + k * n] * B[k + j * n];
          }
          else if((n - 2) % 4 == 0)
          {
              C[(n - 1) + j * n] += A[(n - 1) + k * n] * B[k + j * n];
              C[(n - 2) + j * n] += A[(n - 2) + k * n] * B[k + j * n];
          }
          else if((n - 3) % 4 == 0)
          {
              C[(n - 1) + j * n] += A[(n - 1) + k * n] * B[k + j * n];
              C[(n - 2) + j * n] += A[(n - 2) + k * n] * B[k + j * n];
              C[(n - 3) + j * n] += A[(n - 3) + k * n] * B[k + j * n];
          }
      }
    }
}
