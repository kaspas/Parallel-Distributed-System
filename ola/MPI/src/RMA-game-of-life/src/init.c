#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>
#include <mpi.h>
#include <omp.h>
#include <game-of-life.h>

#define NTHREADS 8

/* set everthing to zero */

void initialize_board (int8_t *board, int N,int M) {
  int   i, j;
  #pragma omp parallel for private(i,j) num_threads(NTHREADS) 
  for (i=0; i<N; i++)
    for (j=0; j<N; j++) 
      Board(i,j) = 0;
}

/* generate random table */

void generate_table (int8_t *board, int N,int M, float threshold) {

  int   i, j;
  //int counter = 0;
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  srand((rank+1)*time(NULL));

  for (i=0; i<N; i++) {

    for (j=0; j<M; j++) {
      Board(i,j) = ( (float)rand() / (float)RAND_MAX ) < threshold;
     // counter += Board(i,j);
    }
  }
}

