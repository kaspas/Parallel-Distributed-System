#include "stdio.h"
#include "stdlib.h"
#include <string.h>
#include "utils.h"

#define MAXBINS 8

int allow_parallelism=0;
int active_threads;
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
  
  if(lv==0){
	  available_threads=numOfThreads;//arxikopoiw ta available threads isa me ta numOfThreads olika threads
	  }
  
  
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
    for(int j=0; j<N; j++){
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
    
	//printf("%d-!",available_threads);
	
    //swap the index pointers  
    swap(&index, &permutation_vector);
    int i;
    //swap the code pointers 
    swap_long(&morton_codes, &sorted_morton_codes);
	if(available_threads>0){//elegxos an yparxoun threads available
	#pragma omp critical//critical section giati ka8orizei an 8a einai nested i oxi dld an ena child borei na dimiourgisei ki alla thread
	{
		//aplos elegxos kai dinw timi sthn allow_parallelism
		if(available_threads>0){
			allow_parallelism=1;
			}
		else{
			allow_parallelism=0;
			}
	}   
	/*if(available_threads!=0){
	printf("%d-!",available_threads);
	}*/
		omp_set_nested(allow_parallelism);
		#pragma omp critical//critical section miwnw ta available threads 
				// den 8a epitrapei sto thread na kanei nested apo to proigoumeno critical section
			  {
			  if(omp_get_nested()!=0 && available_threads>=MAXBINS){
				  active_threads=MAXBINS;
				  available_threads=available_threads-MAXBINS;
				  
				  }
			  else if (omp_get_nested()!=0 && available_threads<MAXBINS && available_threads>=0){
				  active_threads=available_threads;
				  available_threads=0;
				  }
			  /*else{
				  available_threads=0;
				  active_threads=1;
				  }*/
			  }
		#pragma omp parallel private(i) num_threads(active_threads)
		{ //mia parallili for me scheduler static oste na spaei se isa komatia i for
		    #pragma omp for nowait \
		    schedule(static)
			for(i=0; i<MAXBINS; i++){
			    /* if(available_threads!=0){
	   printf("%d-!",available_threads);
	   }*/
			  int offset = (i>0) ? BinSizes[i-1] : 0;
			  int size = BinSizes[i] - offset;
			  
			  truncated_radix_sort(&morton_codes[offset], 
					   &sorted_morton_codes[offset], 
					   &permutation_vector[offset], 
					   &index[offset], &level_record[offset], 
					   size, 
					   population_threshold,
					   sft-3, lv+1);
			/*if(omp_get_nested()!=0 && available_threads>=0){
				  available_threads++;
				  }*/
			
			
			}
		}
	}

	else {//an den yparxoun threads apla trekse siriaka mi peraseis apo ta critical sections

 







    /* Call the function recursively to split the lower levels */
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

