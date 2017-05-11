#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "float.h"
#include "pthread.h"
#include "utils.h"
#define DIM 3

typedef struct{
	unsigned int *codes;
	float *X;
	float *low;
	float step;
	int N;
	int start;	
	int threadID;
	}hash_struct;

inline unsigned int compute_code(float x, float low, float step){

  return floor((x - low) / step);

}


/* Function that does the quantization */
void *thread_quantize_function(void *arg){
  hash_struct *data=(hash_struct *)arg;
  for(int i=0; i<data->N; i++){
    for(int j=0; j<DIM; j++){
      data->codes[data->start+i*DIM + j] = compute_code(data->X[data->start+i*DIM + j], data->low[j], data->step); 
    }
  }
  pthread_exit(NULL);
}
void quantize(unsigned int *codes, float *X, float *low, float step, int N){
	pthread_t *thread;
	hash_struct *arg;
	int flag;
	thread=(pthread_t *)malloc(num_threads*sizeof(pthread_t ));
	arg=(hash_struct *)malloc(num_threads*sizeof(hash_struct));
	for (int i=0;i<num_threads;i++){
		arg[i].threadID=i;//apla gia dokimes
		arg[i].codes =codes;
		arg[i].X=X;
		arg[i].low=low;
		arg[i].step=step;
		arg[i].N=N/num_threads;
		arg[i].start=i*(N/num_threads)*DIM;//starting point
		flag=pthread_create(&thread[i],NULL,thread_quantize_function,(void *)&arg[i]);
		if(flag){
			printf("Error: pthread_create returned code: %d\n", flag);
            return;
            }
		}
	for(int i=0;i<num_threads;i++){
		flag=pthread_join(thread[i],NULL);
		if(flag){
			printf("Error: pthread_join returned code: %d\n", flag);
			return;
		}
	}
    		
}


float max_range(float *x){

  float max = -FLT_MAX;
  for(int i=0; i<DIM; i++){
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

  for(int i=0; i<DIM; i++){
    range[i] = fabs(max[i] - min[i]); // The range of the data
    range[i] += 0.01*range[i]; // Add somthing small to avoid having points exactly at the boundaries 
  }

  qstep = max_range(range) / nbins; // The quantization step 
  
  quantize(codes, X, min, qstep, N); // Function that does the quantization

}



