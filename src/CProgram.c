#include "CInterface.h"

#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>
/*
   Compile and link:

   gcc -c CProgram.c -o CProgram.o
   g++ CProgram.o build/CMakeFiles/DIAlign.dir/*.o -o MyMain_c


*/

int main(int argc, char* argv[]) 
{
  double* r1 = (double*) malloc(12 * sizeof(double));
  double* r2 = (double*) malloc(12 * sizeof(double));
  double* tA = (double*) malloc(2 * sizeof(double));
  double* tB = (double*) malloc(2 * sizeof(double));
  printf("Printing from C!\n");

  memset(r1, 0, 12 * sizeof(double));
  memset(r2, 0, 12 * sizeof(double));

#if 1
  r1[0] = 2;
  r1[1] = 8;
  r2[0] = 9;
  r2[1] = 4;
#endif

  struct CAffineAlignObj res = alignChromatogramsC(6,
                                 r1, 2,
                                 r2, 2,
                                 "none",
                                 tA, 2,
                                 tB, 2,
                                 "mean", "dotProduct",
                                 0.0, 0.0, 0,
                                 0.125, 40,
                                 0.3, true,
                                 0.96, 0.5,
                                 false, 100.0);

  printf("Length of the score: %d \n", res.n_score);
  for (int k = 0; k < res.n_score; k++) 
    printf("Score %d: %f \n", k, res.score[k]);

  printf("Length of index A : %d \n", res.n_indexA_aligned);
  for (int k = 0; k < res.n_indexA_aligned; k++) 
    printf("Index A %d: %d \n", k, res.indexA_aligned[k]);

  printf("Length of index B : %d \n", res.n_indexA_aligned);
  for (int k = 0; k < res.n_indexB_aligned; k++) 
    printf("Index B %d: %d \n", k, res.indexB_aligned[k]);

}

