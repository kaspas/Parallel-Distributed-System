#include "stdio.h"
#include "stdlib.h"
#include "string.h"

#define DIM 3


void data_rearrangement(float *Y, float *X, 
			unsigned int *permutation_vector, 
			int N){
  int i;
  extern int numOfThreads;
  //apla parallili for tis openmp num_threads ka8orizei ton ari8mo twn thread
      #pragma omp parallel num_threads(numOfThreads)
	{
		#pragma omp for private (i) 	

	for(i=0; i<N; i++){
		memcpy(&Y[i*DIM], &X[permutation_vector[i]*DIM], DIM*sizeof(float));
		}
	}
}
