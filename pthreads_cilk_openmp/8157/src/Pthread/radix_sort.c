#include "stdio.h"
#include "stdlib.h"
#include <string.h>
#include "utils.h"
#include "pthread.h"

#define MAXBINS 8

int available_threads;
pthread_mutex_t mut;
typedef struct{
	unsigned long int *morton_codes;
	unsigned long int *sorted_morton_codes; 
	unsigned int *permutation_vector;
	unsigned int *index;
	unsigned int *level_record;
	int N; 
	int population_threshold;
	int sft;
	int lv;
	int offset;
	}radixStruct;
void truncated_radix_sort(unsigned long int *morton_codes, 
			  unsigned long int *sorted_morton_codes, 
			  unsigned int *permutation_vector,
			  unsigned int *index,
			  unsigned int *level_record,
			  int N, 
			  int population_threshold,
			  int sft, int lv);
			  



void first_call_truncated_radix_sort(unsigned long int *morton_codes, 
			  unsigned long int *sorted_morton_codes, 
			  unsigned int *permutation_vector,
			  unsigned int *index,
			  unsigned int *level_record,
			  int N, 
			  int population_threshold,
			  int sft, int lv){
	pthread_mutex_init(&mut,NULL);
	available_threads=num_threads;

	truncated_radix_sort(morton_codes,sorted_morton_codes,permutation_vector,index,level_record,N,population_threshold,sft,lv);
			
}

void *thread_function(void *a){

   
	//printf("%d\n",available_threads);
	// apla vazei orismata stin radix oti exei to radixStruct
	radixStruct *arg=(radixStruct *) a;
	truncated_radix_sort(&arg->morton_codes[arg->offset], 
			   &arg->sorted_morton_codes[arg->offset], 
			   &arg->permutation_vector[arg->offset], 
			   &arg->index[arg->offset], &arg->level_record[arg->offset], 
			   arg->N, 
			   arg->population_threshold,
			   arg->sft, arg->lv);
		/*pthread_mutex_lock(&mut);
		available_threads++;
		pthread_mutex_unlock(&mut);*/
		pthread_exit(NULL);
	}
inline void swap_long(unsigned long int **x, unsigned long int **y){

  unsigned long int *tmp;
  tmp = x[0];
  x[0] = y[0];
  y[0] = tmp;

}

inline void swap(unsigned int **x, unsigned int **y){

  unsigned int *tmp;
  tmp = x[0];
  x[0] = y[0];
  y[0] = tmp;

}

void truncated_radix_sort(unsigned long int *morton_codes, 
			  unsigned long int *sorted_morton_codes, 
			  unsigned int *permutation_vector,
			  unsigned int *index,
			  unsigned int *level_record,
			  int N, 
			  int population_threshold,
			  int sft, int lv){
  int BinSizes[MAXBINS] = {0};
  int BinCursor[MAXBINS] = {0};
  unsigned int *tmp_ptr;
  unsigned long int *tmp_code;
  /*if(lv==0){
	pthread_mutex_init(&mut,NULL);
	available_threads=num_threads;
	  }*/
  //if(lv>=10)printf("%d\n",lv);
  if(N<=0){
    return;
  }
  else if(N<=population_threshold || sft < 0) { // Base case. The node is a leaf
    level_record[0] = lv; // record the level of the node
    memcpy(permutation_vector, index, N*sizeof(unsigned int)); // Copy the pernutation vector
    memcpy(sorted_morton_codes, morton_codes, N*sizeof(unsigned long int)); // Copy the Morton codes 
    return;
  }
  else{
    level_record[0] = lv;
    // Find which child each point belongs to 
    for(int j=0; j<N;
     j++){
      unsigned int ii = (morton_codes[j]>>sft) & 0x07;
      BinSizes[ii]++;
    }
    // scan prefix (must change this code)  
    int offset = 0;
    for(int i=0; i<MAXBINS; i++){
      int ss = BinSizes[i];
      BinCursor[i] = offset;
      offset += ss;
      BinSizes[i] = offset;
    }
    
    for(int j=0; j<N; j++){
      unsigned int ii = (morton_codes[j]>>sft) & 0x07;
      permutation_vector[BinCursor[ii]] = index[j];
      sorted_morton_codes[BinCursor[ii]] = morton_codes[j];
      BinCursor[ii]++;
    }
    
    //swap the index pointers  
    swap(&index, &permutation_vector);

    //swap the code pointers 
    swap_long(&morton_codes, &sorted_morton_codes);

    /* Call the function recursively to split the lower levels */
    int flag;
    if(available_threads>0){
		//paralilopoihsh tis anadromis einai critical na min boun alla thread edw
		//opote xrisimopoiw mutex
		pthread_mutex_lock(&mut);
		if(available_threads>0 && available_threads>=MAXBINS){//to available_threads>0 perito ksaxasa na to svisw
			pthread_t *thread=(pthread_t *)malloc(MAXBINS*sizeof(pthread_t ));
			radixStruct *arg=(radixStruct *)malloc(MAXBINS*sizeof(radixStruct ));
			available_threads=available_threads-MAXBINS;
			pthread_mutex_unlock(&mut);
			for(int i=0; i<MAXBINS; i++){
			  int offset = (i>0) ? BinSizes[i-1] : 0;
			  int size = BinSizes[i] - offset;
			  //available_threads--;
			  arg[i].morton_codes=morton_codes;
			  arg[i].sorted_morton_codes=sorted_morton_codes;
			  arg[i].permutation_vector=permutation_vector;
			  arg[i].index=index;
			  arg[i].level_record=level_record;
			  arg[i].N=size;
			  arg[i].population_threshold=population_threshold;
			  arg[i].sft=sft-3;
			  arg[i].lv=lv+1;
			  arg[i].offset=offset;
			  flag=pthread_create(&thread[i],NULL,thread_function,(void *)(&arg[i]));
			 //printf("Bika edw!!\n");
			if(flag){//elegxos an egine create
			  printf("Error: pthread_create returned code: %d\n", flag);
			  return;
				}

			}
			  //pthread_mutex_unlock(&mut);

			for(int i=0; i<MAXBINS; i++){
				flag=pthread_join(thread[i],NULL);
			if(flag){//elegxos an egine join
			  printf("Error: pthread_join returned code: %d\n", flag);
			  return;				
				}
			}
			free(arg);
					free(thread);

		}
		else if(available_threads>0 && available_threads<MAXBINS) {//to available_threads>0 perito ksaxasa na to svisw
			int reps=available_threads;//apo8ikeuw to available threads giati 8a to midenisw
			//printf("Bika edw!!\n");
			pthread_t *thread=(pthread_t *)malloc(available_threads*sizeof(pthread_t ));
			radixStruct *arg=(radixStruct *)malloc(available_threads*sizeof(radixStruct ));
		    available_threads=0;
		    pthread_mutex_unlock(&mut);
			for(int i=0; i<reps; i++){
			  int offset = (i>0) ? BinSizes[i-1] : 0;
			  int size = BinSizes[i] - offset;
			  //available_threads--;
			  arg[i].morton_codes=morton_codes;
			  arg[i].sorted_morton_codes=sorted_morton_codes;
			  arg[i].permutation_vector=permutation_vector;
			  arg[i].index=index;
			  arg[i].level_record=level_record;
			  arg[i].N=size;
			  arg[i].population_threshold=population_threshold;
			  arg[i].sft=sft-3;
			  arg[i].lv=lv+1;
			  arg[i].offset=offset;
			  flag=pthread_create(&thread[i],NULL,thread_function,(void *)(&arg[i]));
			  if(flag){//elegxos an egine create
				printf("Error: pthread_create returned code: %d\n", flag);
				}
			}


			//pthread_mutex_unlock(&mut);	

			//sinexizw tis ypoloipes anadromes siriaka
			for(int i=reps; i<MAXBINS; i++){
				int offset = (i>0) ? BinSizes[i-1] : 0;
				int size = BinSizes[i] - offset;
		  
				truncated_radix_sort(&morton_codes[offset], 
				   &sorted_morton_codes[offset], 
				   &permutation_vector[offset], 
				   &index[offset], &level_record[offset], 
				   size, 
				   population_threshold,
				   sft-3, lv+1);
			}
		   for(int i=0; i<reps; i++){
				flag=pthread_join(thread[i],NULL);
				//printf("Bika edw!!\n");
			if(flag){//elegxos an egine join
			  printf("Error: pthread_join returned code: %d\n", flag);
				}
			}
			free(arg);
			free(thread);
		}
		else if(available_threads<=0){// periptwsi na bei katala8os thread edw enw den yparxoun available threads
						// logo twn mutex prp na to provlepsw
			//	printf("Bika\n");
		pthread_mutex_unlock(&mut);	
		//sinexizei seiriaka
		for(int i=0; i<MAXBINS; i++){
				int offset = (i>0) ? BinSizes[i-1] : 0;
				int size = BinSizes[i] - offset;
		  
				truncated_radix_sort(&morton_codes[offset], 
				   &sorted_morton_codes[offset], 
				   &permutation_vector[offset], 
				   &index[offset], &level_record[offset], 
				   size, 
				   population_threshold,
				   sft-3, lv+1);
			}
		}
	 }
	  else{//an den exw available threads pane seiriaka apefige ta mutex diladi
		for(int i=0; i<MAXBINS; i++){
				int offset = (i>0) ? BinSizes[i-1] : 0;
				int size = BinSizes[i] - offset;
		  
				truncated_radix_sort(&morton_codes[offset], 
				   &sorted_morton_codes[offset], 
				   &permutation_vector[offset], 
				   &index[offset], &level_record[offset], 
				   size, 
				   population_threshold,
				   sft-3, lv+1);
			}		  
		  }
  } 
}

