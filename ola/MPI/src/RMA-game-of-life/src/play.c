#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>
#include <omp.h>
#include <mpi.h>
#include <game-of-life.h>
#define NTHREADS 8
void play (int8_t *board, int8_t *newboard, int N, int M) {
  /*
    (copied this from some web page, hence the English spellings...)

    1.STASIS : If, for a given cell, the number of on neighbours is 
    exactly two, the cell maintains its status quo into the next 
    generation. If the cell is on, it stays on, if it is off, it stays off.

    2.GROWTH : If the number of on neighbours is exactly three, the cell 
    will be on in the next generation. This is regardless of the cell's
    current state.

    3.DEATH : If the number of on neighbours is 0, 1, 4-8, the cell will 
    be off in the next generation.
  */

  int rank,nproc;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);
  int   i, j, a;

  if(nproc==1){

	  /* for each cell, apply the rules of Life */
	#pragma omp parallel for schedule(dynamic) private(i,j,a) shared(N,M,board,newboard) num_threads(NTHREADS) 
	  for (i=0; i<N; i++)
		for (j=0; j<M; j++) {
		  a = adjacent_to (board, i, j, N,M);
		  if (a == 2) NewBoard(i,j) = Board(i,j);
		  if (a == 3) NewBoard(i,j) = 1;
		  if (a < 2) NewBoard(i,j) = 0;
		  if (a > 3) NewBoard(i,j) = 0;
		}

	  /* copy the new board back into the old board */

	}
  else{
	  int8_t first_row[M],last_row[M],local_fr[M],local_lr[M];
	  MPI_Win frwin,lrwin;
	  MPI_Win_create(first_row,M,sizeof(int8_t),MPI_INFO_NULL,MPI_COMM_WORLD,&frwin);
	  MPI_Win_create(last_row,M,sizeof(int8_t),MPI_INFO_NULL,MPI_COMM_WORLD,&lrwin);
	  MPI_Alloc_mem(M*sizeof(int8_t),MPI_INFO_NULL,first_row);
	  MPI_Alloc_mem(M*sizeof(int8_t),MPI_INFO_NULL,last_row);
	  #pragma omp parallel for schedule(dynamic) private (i) shared(board)
      for(i=0;i<M;i++){
		first_row[i]=Board(0,i);
		last_row[i]=Board(N-1,i);
      }

	  MPI_Win_fence(0,frwin);
	  MPI_Win_fence(0,lrwin);
	  
      int next=rank+1;
      int prev=rank-1;
      
      if(rank==0) prev=nproc-1;
      if(rank==nproc-1) next=0;
	  /* for each cell, apply the rules of Life */
	 #pragma omp parallel for schedule(dynamic) private(i,j,a) shared(newboard,board,M) num_threads(NTHREADS) 
	  for (i=1; i<N-1; i++)
		for (j=0; j<M; j++) {
		  a = adjacent_to (board, i, j, N-1, M);
		  if (a == 2) NewBoard(i,j) = Board(i,j);
		  if (a == 3) NewBoard(i,j) = 1;
		  if (a < 2) NewBoard(i,j) = 0;
		  if (a > 3) NewBoard(i,j) = 0;
		}

	  /* copy the new board back into the old board */


	MPI_Get(&local_fr[0],M,MPI_INT8_T,prev,0,M,MPI_INT8_T,frwin);
    MPI_Get(&local_lr[0],M,MPI_INT8_T,next,0,M,MPI_INT8_T,lrwin);	
	MPI_Win_fence(0,frwin);
	MPI_Win_fence(0,lrwin);
	/*if(rank==0){
		for(int i=0;i<M;i++){
			printf("%d-",local_fr[i]);
			}
		printf("\n");
		}*/
	int8_t *temp1,*temp2;
	temp1=(int8_t *)malloc(M*3*sizeof(int8_t));
	temp2=(int8_t *)malloc(M*3*sizeof(int8_t));
	#pragma omp parallel for schedule(dynamic) private(j) shared(board,local_fr,local_lr,temp1,temp2,M) num_threads(NTHREADS) 
	for(j=0;j<M;j++){
		temp1[j]=local_fr[j];
		temp1[M+j]=Board(0,j);
		temp1[2*M+j]=Board(1,j);
		temp2[j]=Board(N-2,j);
		temp2[M+j]=Board(N-1,j);
		temp2[2*M+j]=local_lr[j];
		}
	#pragma omp parallel for schedule(dynamic) private(j,a) shared(newboard,board,M) num_threads(NTHREADS) 	
	for(j=0;j<M;j++){
		a = adjacent_to(temp1,1,j,3, M);
		if (a == 2) NewBoard(0,j) = Board(0,j);
		if (a == 3) NewBoard(0,j) = 1;
		if (a < 2) NewBoard(0,j) = 0;
		if (a > 3) NewBoard(0,j) = 0;
		}
	#pragma omp parallel for schedule(dynamic) private(j,a) shared(newboard,board,M) num_threads(NTHREADS) 
	for (j=0; j<M; j++) {
		  a = adjacent_to (temp2, 1, j, 3, M);
		  if (a == 2) NewBoard(N-1,j) = Board(N-1,j);
		  if (a == 3) NewBoard(N-1,j) = 1;
		  if (a < 2) NewBoard(N-1,j) = 0;
		  if (a > 3) NewBoard(N-1,j) = 0;
	  }
	free(temp1);
	free(temp2);
	MPI_Win_free(&frwin);
	MPI_Win_free(&lrwin);
  }
}
