/*
 * Game of Life implementation based on
 * http://www.cs.utexas.edu/users/djimenez/utsa/cs1713-3/c/life.txt
 * 
 */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>
#include <omp.h>
#include <mpi.h>
#include <game-of-life.h>
#include <sys/time.h>
#define NTHREADS 8

int main (int argc, char *argv[]) {
  int8_t   *board, *newboard; 
  int i,rc;
  struct timeval startwtime, endwtime;
  if (argc != 6) { // Check if the command line arguments are correct 
    printf("Usage: %s N thres disp\n"
	   "where\n"
	   "  N     : size of table (N x N)\n"
	   "  thres : propability of alive cell\n"
           "  t     : number of generations\n"
	   "  disp  : {1: display output, 0: hide output}\n"
	   "  M 	: Number of Columns \n	"
           , argv[0]);
    return (1);
  }
  int nproc,rank;
  int len;
  char name[MPI_MAX_PROCESSOR_NAME];
  rc=MPI_Init(&argc,&argv);
  if (rc != MPI_SUCCESS) {
    printf("Error starting MPI program.\nTerminating.\n");
    MPI_Abort(MPI_COMM_WORLD, rc);
  } 
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  MPI_Comm_size(MPI_COMM_WORLD,&nproc);
  MPI_Get_processor_name(name,&len);
  printf("Hi from process :%d from processor : %s\n",rank,name);
  omp_set_num_threads(NTHREADS);
  //Input command line arguments
  int M=atoi(argv[5]);
  int N = atoi(argv[1]);        // Array size
  double thres = atof(argv[2]); // Propability of life cell
  int t = atoi(argv[3]);        // Number of generations 
  int disp = atoi(argv[4]);     // Display output?
  N=N/nproc;
  MPI_Barrier(MPI_COMM_WORLD);
 if(rank==0){
	  gettimeofday(&startwtime,NULL);
	  }
  printf("Size %dx%d with propability: %0.1f%%\n", N, M, thres*100);
  
  board = NULL;
  newboard = NULL;
  board = (int8_t *)malloc(N*M*sizeof(int8_t));
  
  if (board == NULL){
    printf("\nERROR: Memory allocation did not complete successfully!\n");
    MPI_Finalize();
    return (1);
  }

  /* second pointer for updated result */
  newboard = (int8_t *)malloc(N*M*sizeof(int8_t));

  if (newboard == NULL){
    printf("\nERROR: Memory allocation did not complete successfully!\n");
    MPI_Finalize();
    return (1);
  }
  //int j;
  //initialize_board (board, N,M);
  printf("Board initialized\n");
  generate_table (board, N, M, thres);
  printf("Board generated\n");
  if(rank==0){
	gettimeofday(&endwtime,NULL);
	double genTime=(double)((endwtime.tv_usec - startwtime.tv_usec)
			/1.0e6 + endwtime.tv_sec - startwtime.tv_sec);
			printf("%lf\n",genTime);
	}

  /* play game of life 100 times */

  for (i=0; i<t; i++) {
    if (disp) display_table (board, N, M);
    play (board, newboard, N, M);  
    int8_t *temp=board;	  
	board=newboard;  
    newboard=temp;
  }
    printf("Game finished after %d generations.\n", t);

   MPI_Barrier(MPI_COMM_WORLD);
   if(rank==0){
	 gettimeofday(&endwtime,NULL);
	 double executionTime=(double)((endwtime.tv_usec - startwtime.tv_usec)
				/1.0e6 + endwtime.tv_sec - startwtime.tv_sec);
				printf("%lf\n",executionTime);
	  }
   free(board);
   free(newboard);
   MPI_Finalize();
   return(0);  
}


