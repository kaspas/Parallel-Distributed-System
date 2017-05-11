#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "utils.h"
#include "pthread.h"
#define DIM 3

typedef struct {
	unsigned long int *mcodes;
	unsigned int *codes;
	int N;
	int start;
	int threadID;
	}morton_struct;


inline unsigned long int splitBy3(unsigned int a){
    unsigned long int x = a & 0x1fffff; // we only look at the first 21 bits
    x = (x | x << 32) & 0x1f00000000ffff;  // shift left 32 bits, OR with self, and 00011111000000000000000000000000000000001111111111111111
    x = (x | x << 16) & 0x1f0000ff0000ff;  // shift left 32 bits, OR with self, and 00011111000000000000000011111111000000000000000011111111
    x = (x | x << 8) & 0x100f00f00f00f00f; // shift left 32 bits, OR with self, and 0001000000001111000000001111000000001111000000001111000000000000
    x = (x | x << 4) & 0x10c30c30c30c30c3; // shift left 32 bits, OR with self, and 0001000011000011000011000011000011000011000011000011000100000000
    x = (x | x << 2) & 0x1249249249249249;
    return x;
}

inline unsigned long int mortonEncode_magicbits(unsigned int x, unsigned int y, unsigned int z){
    unsigned long int answer;
    answer = splitBy3(x) | splitBy3(y) << 1 | splitBy3(z) << 2;
    return answer;
}

void *thread_morton_function(void *a){
	morton_struct *arg=(morton_struct *)a;
	
	for(int i=0; i<arg->N; i++){
    // Compute the morton codes from the hash codes using the magicbits mathod
    arg->mcodes[arg->start+i] = mortonEncode_magicbits(arg->codes[(arg->start)*DIM+i*DIM], arg->codes[(arg->start)*DIM+i*DIM + 1], arg->codes[(arg->start)*DIM+i*DIM + 2]);
   }
   pthread_exit(NULL);
}

/* The function that transform the morton codes into hash codes */ 
void morton_encoding(unsigned long int *mcodes, unsigned int *codes, int N, int max_level){
    int flag;
    pthread_t *thread;
    thread=(pthread_t *)malloc(num_threads*sizeof(pthread_t ));  
    morton_struct *arg;
    arg=(morton_struct *)malloc(num_threads*sizeof(morton_struct));
    for(int i=0;i<num_threads;i++){
		arg[i].threadID=i;//apla gia dokimes
		arg[i].codes=codes;
		arg[i].mcodes=mcodes;
		arg[i].N=N/num_threads;
		arg[i].start=i*(N/num_threads);//starting point
		flag=pthread_create(&thread[i],NULL,thread_morton_function,(void *)&arg[i]);
		if (flag){
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


