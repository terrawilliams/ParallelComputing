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

#define BLOCK_SIZE 50

#define min(a,b) (((a)<(b))?(a):(b))

const char* dgemm_desc = "Naive, three-loop dgemm.";

/* This routine performs a dgemm operation
 *  C := C + A * B
 * where A, B, and C are lda-by-lda matrices stored in column-major format.
 * On exit, A and B maintain their input values. */    
void square_dgemm (int n, double* A, double* B, double* restrict C)
{
    double AT[n * n];
//#pragma loop_count min=31, max=769, avg=345
    for(int i1 = 0; i1 < n; i1 += BLOCK_SIZE)
    {
        for(int j1 = 0; j1 < n; j1 += BLOCK_SIZE)
        {
            int i_edge = min(BLOCK_SIZE, n - i);
            double * AT2 = AT + j1 * n + i1;
            for (int i = 0; i < i_edge; i++)
            {
#pragma ivdep
                for (int j = j1; j < min(n, j1 + BLOCK_SIZE); j++)
                {
                    AT2[j * n + i] = A[i * n + j];
                }
            }
        }
    }

        for (int j = 0; j < n; ++j)
        {
            for (int i = 0; i < n; i++)
            {
#pragma vector unaligned
                for (int k = 0; k < n; k++)
                {
                    C[i + j * n] += AT[k + i * n] * B[k + j * n];
                }
            }
        }
}
