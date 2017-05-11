#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "utils.h"
#include "pthread.h"

#define DIM 3

typedef struct{
	float *Y;
	float *X;
	unsigned int *permutation_vector;
	int start;
	int N;
	int threadID;
	}data_rearrangement_struct;
	
	
void *thread_rearrangement(void *a){
  data_rearrangement_struct *arg=(data_rearrangement_struct *)a;
  for(int i=0; i<arg->N; i++){
    memcpy(&arg->Y[(arg->start)*DIM+i*DIM], &arg->X[arg->permutation_vector[arg->start+i]*DIM], DIM*sizeof(float));
  }	
}

void data_rearrangement(float *Y, float *X, 
			unsigned int *permutation_vector, 
			int N){
    pthread_t *thread;
    data_rearrangement_struct *arg;
    int flag;
	thread=(pthread_t *)malloc(num_threads*sizeof(pthread_t ));
	arg=(data_rearrangement_struct *)malloc(num_threads*sizeof(data_rearrangement_struct ));
	
    for (int i=0;i<num_threads;i++){
	    arg[i].Y=Y;
	    arg[i].X=X;
	    arg[i].permutation_vector=permutation_vector;
	    arg[i].start=i*(N/num_threads);//starting point
	    arg[i].threadID;//apla gia dokimes
	    arg[i].N=N/num_threads;
	    flag=pthread_create(&thread[i],NULL,thread_rearrangement,(void *)&arg[i]);
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



