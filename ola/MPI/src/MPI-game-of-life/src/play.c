#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>
#include <game-of-life.h>
#include <mpi.h>
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
  int i,j, a;
  int next,prev,nproc,rank;
  int8_t *recv_buffer1,*recv_buffer2,*send_buffer1,*send_buffer2;
  recv_buffer1=(int8_t *)malloc(M*sizeof(int8_t));
  recv_buffer2=(int8_t *)malloc(M*sizeof(int8_t));
  send_buffer1=(int8_t *)malloc(M*sizeof(int8_t));
  send_buffer2=(int8_t *)malloc(M*sizeof(int8_t));
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);
  send_buffer1=&Board(0,0);
  send_buffer2=&Board(N-1,0);
  MPI_Request send_req[2],recv_req[2];  
  prev=rank-1;
  next=rank+1;
  //MPI Communication 
  if(rank==0) prev=nproc-1;
  if(rank==(nproc-1)) next=0;
  if(nproc==1){

			/* for each cell, apply the rules of Life */
	  #pragma omp parallel for schedule(dynamic) private(i,j,a) shared(N,M,board,newboard) collapse(2) num_threads(NTHREADS) 
	  for (i=0; i<N; i++){
		  for (j=0; j<M; j++) {
		  a = adjacent_to (board, i, j, N, M);
		  if (a == 2) NewBoard(i,j) = Board(i,j);
		  if (a == 3) NewBoard(i,j) = 1;
		  if (a < 2) NewBoard(i,j) = 0;
		  if (a > 3) NewBoard(i,j) = 0;
	  }
	}
	  /* copy the new board back into the old board */
	  
	//  #pragma omp parallel for schedule(dynamic) private(i,j) shared(newboard,board) num_threads(NTHREADS) 

//	  for(i=0;i<N;i++){
	//	for(j=0;j<M;j++)
	//	Board(i,j)=NewBoard(i,j);
		  
	   // }
  }    
  else{
  	  MPI_Irecv(recv_buffer1,M,MPI_INT8_T,prev,prev,MPI_COMM_WORLD,&recv_req[0]);
	  MPI_Irecv(recv_buffer2,M,MPI_INT8_T,next,next,MPI_COMM_WORLD,&recv_req[1]);
 	  MPI_Isend(send_buffer1,M,MPI_INT8_T,prev,rank,MPI_COMM_WORLD,&send_req[0]);
	  MPI_Isend(send_buffer2,M,MPI_INT8_T,next,rank,MPI_COMM_WORLD,&send_req[1]);

	  /* for each cell, apply the rules of Life */
	  #pragma omp parallel for schedule(dynamic) private(i,j,a) shared(newboard,board,N,M) num_threads(NTHREADS) 
	  for (i=1; i<N-1; i++){
		  for (j=0; j<M; j++) {
		  a = adjacent_to (board, i, j, N-1, M);
		  if (a == 2) NewBoard(i,j) = Board(i,j);
		  if (a == 3) NewBoard(i,j) = 1;
		  if (a < 2) NewBoard(i,j) = 0;
		  if (a > 3) NewBoard(i,j) = 0;
	  }
	}
	  /* copy the new board back into the old board */
	  
	/*#pragma omp parallel for schedule(dynamic) private(i,j) shared(newboard,board,N,M) num_threads(NTHREADS) 
	  for (i=1; i<N-1; i++){
		for (j=0; j<M; j++) 
		  Board(i,j)=NewBoard(i,j);
	  }*/
	  MPI_Status status[2];
	  MPI_Waitall(2,send_req,status	);
	  MPI_Waitall(2,recv_req,status);
	  
	  int8_t *temp,*temp2;
	  temp=(int8_t *)malloc(3*M*sizeof(int8_t ));
	  temp2=(int8_t *)malloc(3*M*sizeof(int8_t ));
	 #pragma omp parallel for schedule(dynamic) private(i) shared(board,recv_buffer1,recv_buffer2,temp,temp2,M) num_threads(NTHREADS) 
	  for(i=0;i<M;i++){
		  temp[i]=recv_buffer1[i];
		  temp[M+i]=Board(0,i);
		  temp[2*M+i]=Board(1,i);
		  temp2[i]=Board(N-2,i);
		  temp2[1*M+i]=Board(N-1,i);
		  temp2[2*M+i]=recv_buffer2[i];
		  }
	#pragma omp parallel for schedule(dynamic) private(j,a) shared(newboard,board,M) num_threads(NTHREADS) 
      for (j=0; j<M; j++) {
		  a = adjacent_to (temp, 1, j, 3, M);
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
	 /*#pragma omp parallel for schedule(dynamic) private(j) shared(newboard,board,M) num_threads(NTHREADS) 
		for (j=0; j<M; j++) {
		  Board(0,j)=NewBoard(0,j);
		  Board(N-1,j)=NewBoard(N-1,j);
		}*/
	free(temp2);
	free(temp);
    
  }
}
