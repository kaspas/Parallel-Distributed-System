#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "float.h"
#include "utils.h"
#include <cilk/cilk.h>

#define DIM 3

inline unsigned int compute_code(float x, float low, float step){

  return floor((x - low) / step);

}


/* Function that does the quantization */
/*void quantize(unsigned int *codes, float *X, float *low, float step, int N){
    for(int j=0; j<DIM; j++){
    cilk_for(int i=0; i<N; i++){
      codes[i*DIM + j] = compute_code(X[i*DIM + j], low[j], step); 
    }
  }

}*/

void quantize(unsigned int *codes, float *X, float *low, float step, int N){
	
         int i;
    cilk_for(i=0; i<N; i++){
      codes[i*DIM ] = compute_code(X[i*DIM ], low[0], step); 
      codes[i*DIM + 1] = compute_code(X[i*DIM + 1], low[1], step); 
      codes[i*DIM + 2] = compute_code(X[i*DIM + 2], low[2], step); 
    }
}

float max_range(float *x){

  float max = -FLT_MAX;
  int i;
  for(i=0; i<DIM; i++){
    if(max<x[i]){
      max = x[i];
    }
  }

  return max;

}

void compute_hash_codes(unsigned int *codes, float *X, int N, 
			int nbins, float *min, 
			float *max){
  
  float range[DIM];
  float qstep;
  int i;
  for(i=0; i<DIM; i++){
    range[i] = fabs(max[i] - min[i]); // The range of the data
    range[i] += 0.01*range[i]; // Add somthing small to avoid having points exactly at the boundaries 
  }

  qstep = max_range(range) / nbins; // The quantization step 

  quantize(codes, X, min, qstep, N); // Function that does the quantization
}



